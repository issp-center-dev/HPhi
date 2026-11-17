"""Explicit unavailable/unknown backends must fail before running HPhi or MPI."""
from pathlib import Path
import os
import subprocess
import sys
import tempfile


root = Path(tempfile.mkdtemp(prefix='fd-backend-requirements-', dir=Path.cwd()))
scripts = ('symmetry_nconserved_core.py', 'fulldiag_spinless_backends.py')
cases = (
    ('missing_scalapack', 'ScaLAPACK', 'missing-mpiexec -np 2', {}, 'required backends unavailable'),
    ('missing_elpa', 'LAPACK,ScaLAPACK,ELPA', 'missing-mpiexec -np 2',
     {'HPHI_HAS_SCALAPACK': '1'}, 'required backends unavailable'),
    ('single_rank', 'ScaLAPACK,ELPA', 'missing-mpiexec -np 1',
     {'HPHI_HAS_SCALAPACK': '1', 'HPHI_HAS_ELPA': '1'}, 'required backends unavailable'),
    ('unknown', 'ELAP', '', {}, 'unknown required backends'),
)
for script in scripts:
    for name, required, launcher, flags, diagnostic in cases:
        path = root/f'{Path(script).stem}_{name}'
        path.mkdir()
        env = {k: v for k, v in os.environ.items()
               if k not in ('HPHI_HAS_SCALAPACK', 'HPHI_HAS_ELPA')}
        env.update(flags, HPHI_FD_REQUIRE_BACKENDS=required, MPIRUN=launcher,
                   HPHI_TEST_NCONSERVED_CORE_SCOPE='fulldiag')
        result = subprocess.run([sys.executable, str(Path(__file__).with_name(script)),
                                 'missing-HPhi', 'missing-probe'], cwd=path, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, timeout=30)
        (path/'run.log').write_text(result.stdout)
        assert result.returncode != 0 and diagnostic in result.stdout, (name, result.stdout)
        assert 'FileNotFoundError' not in result.stdout, result.stdout
        assert not (path/'symmetry_nconserved_core').exists()
        assert not list(path.glob('sector-fd-spinless-backends-*'))
        print(f'{Path(script).stem}: {name} rejected before HPhi/MPI', flush=True)
