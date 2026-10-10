"""Collective acceptance on literal 6/14/4-dimensional Kondo sectors.

Dropping empty-rank participation, truncating payloads, using raw moments,
or applying the wrong projected Hamiltonian must fail these comparisons.
"""
import os
from pathlib import Path
import shlex
import shutil
import tempfile
import numpy as np
from symmetry_kondo_common import (expert_empty_rank_case, selected_layouts,
                                   read_manifest, run_probe, read_vector_parts,
                                   read_kondo_checkpoint)
from symmetry_kondo_reference import make_empty_rank_reference, sector_hamiltonian
from symmetry_spingc_common import mpi_size, run_case, checkpoint_digest


def check_ownership(path, stem, dimension, ranks):
    files = sorted(path.glob(stem+'*.dat'))
    assert len(files) == ranks
    empty = 0
    for file in files:
        info = read_manifest(file.with_suffix('.info'))
        size = int(info['local_dim'])
        assert int(info['owned_basis_elements']) == size
        assert int(info['dim']) == dimension
        assert len(file.read_text().splitlines()) == size
        if size == 0:
            empty += 1
            # No owned rows, but the index-zero sentinel and replicated
            # global metadata remain allowed by the storage contract.
            assert not file.read_text()
    assert empty == max(0, ranks-dimension), (empty, ranks, dimension)


def checkpoint_vector(path, dimension, ranks, label='zvo_eigenvec_0'):
    files = [path/'output'/f'{label}_rank_{r}.dat' for r in range(ranks)]
    parts = [read_kondo_checkpoint(f) for f in files]
    assert all(tuple(h[26:28]) == checkpoint_digest(parts) for h, _ in parts)
    assert sum(len(v) == 0 for _, v in parts) == max(0, ranks-dimension)
    for rank, (header, payload) in enumerate(parts):
        assert tuple(header[11:13]) == (ranks, rank)
        assert len(payload) == header[15]
    vector = np.concatenate([v for _, v in parts])
    assert len(vector) == dimension and abs(np.vdot(vector, vector)-1) < 1e-8
    return vector


def requests_and_values(reference, vector):
    # E00 at a local site is a projector, so every positive power equals
    # E00. This checks all output formats without a cell-chain assumption.
    projector = np.array([bool(int(w)&1) for w in reference.words])
    value = np.vdot(reference.basis@vector, projector*(reference.basis@vector))
    row = [0, 0, 0, 0]
    requests = {key: [row*n] for key, n in [('one', 1), ('two', 2),
                ('three', 3), ('four', 4), ('six', 6)]}
    requests['nbody'] = [[5]+row*5]
    return requests, {key: np.array([value]) for key in requests}


def check_correlations(path, reference, vector, suffix):
    from symmetry_kondo_observables import CORRELATION_FILES
    from symmetry_kondo_thermal import finite_data, check_correlation_block
    requests, values = requests_and_values(reference, vector)
    for kind, rows in requests.items():
        file = path/'output'/f'zvo_{CORRELATION_FILES[kind]}{suffix}.dat'
        check_correlation_block(finite_data(file), rows, values[kind], 0, str(file))


def check_manifests(paths, model, layout, dimension, ranks, omp_threads):
    expected = [path/'output/zvo_symmetry_sector.dat' for path in paths]
    assert len(set(expected)) == len(expected)
    for manifest in expected:
        assert manifest.is_file(), ('missing symmetry manifest', manifest)
        metadata = read_manifest(manifest)
        assert metadata['model'] == model
        assert metadata['basis_layout'] == layout
        assert int(metadata['sector_dim']) == dimension
        assert int(metadata['mpi_ranks']) == ranks
        if omp_threads is not None:
            assert int(metadata['omp_threads']) == int(omp_threads)


def run_empty_rank_suite(kind: str, executable: Path, probe: Path | None = None) -> None:
    from symmetry_kondo_observables import (add_correlation_requests, parse_moments,
                                            check_eigenvector)
    from symmetry_kondo_thermal import (captured_initials, tpq_step, ctpq_step,
                                        check_data, NUM_AVE, NSTEPS)
    from symmetry_kondo_te import attach_seed, te_step, exact_endpoint_step
    from symmetry_spingc_te import schedule, TIMES
    launcher = shlex.split(os.environ.get('MPIRUN', ''))
    ranks = mpi_size(launcher)
    assert ranks in (1, 2, 4, 16)
    probe = probe or executable.parent/'unittest_symmetry_sector_probe'
    root = Path(tempfile.mkdtemp(prefix=f'symmetry_kondo_empty_{kind}_', dir='.'))
    print('artifacts:', root.resolve(), flush=True)
    errors = dict(matrix=0., residual=0., moments=0., vector=0., data=0., endpoint=0.)
    for layout in selected_layouts():
        env = dict(os.environ, HPHI_SYMMETRY_BASIS_LAYOUT=layout,
                   HPHI_TEST_SYMMETRY_CAPTURE='1')
        for model, dimension in [('Kondo', 6), ('KondoNConserved', 14), ('KondoGC', 4)]:
            executed_case_paths = []
            ref = make_empty_rank_reference(model)
            assert ref.basis.shape[1] == dimension
            np.testing.assert_allclose(ref.basis.conj().T@ref.basis, np.eye(dimension), atol=1e-12, rtol=0)
            h = sector_hamiltonian(ref)
            assert np.linalg.norm(h-h.conj().T) < 1e-12*max(1, np.linalg.norm(h))
            hb = np.column_stack([ref.apply_h(v) for v in ref.basis.T])
            assert np.linalg.norm(hb-ref.basis@h) < 1e-12*max(1, np.linalg.norm(h))
            prefix = root/f'{model}_{layout}'
            vector = np.arange(1, dimension+1)+1j*np.arange(dimension, 0, -1)
            vector /= np.linalg.norm(vector)
            requests, _ = requests_and_values(ref, vector)
            if kind == 'basis':
                path = expert_empty_rank_case(prefix, model, method=3, options={})
                executed_case_paths.append(path)
                actual = run_probe(path, None, probe, action='matvec', layout=layout)
                np.testing.assert_array_equal(actual['representatives'], ref.representatives)
                error = np.linalg.norm(actual['matrix']-h)
                errors['matrix'] = max(errors['matrix'], error)
                assert error < 1e-11*max(1, np.linalg.norm(h))
                empty = [info for info in actual['info'] if int(info['local_dim']) == 0]
                assert len(empty) == max(0, ranks-dimension)
                assert all(int(info.get('owned_basis_elements', -1)) == int(info['local_dim']) for info in actual['info']), 'missing owned basis evidence'
            elif kind == 'observables':
                path = expert_empty_rank_case(prefix, model, method=3,
                       options={'CalcMod': {'OutputEigenVec': 1}, 'ModPara': {'LanczosEps': 18}})
                executed_case_paths.append(path)
                result = run_probe(path, None, probe, action='moments', layout=layout, vector=vector)
                values = parse_moments(result['text'])
                for key in ('N', 'N2', 'D', 'D2', 'Sz', 'Sz2'):
                    error = abs(values[key]-abs(ref.basis@vector)**2@ref.moments[key])
                    errors['moments'] = max(errors['moments'], error)
                    assert error < 1e-8
                add_correlation_requests(path, requests, False)
                run_case(path, probe, 'cg', launcher, env)
                actual = read_vector_parts(path, 'sector_final_sample0_step0.rank', dimension)
                energy = check_eigenvector(actual, h, np.linalg.eigvalsh(h)[0])
                errors['residual'] = max(errors['residual'], np.linalg.norm(h@actual-energy*actual))
                check_ownership(path, 'sector_final_sample0_step0.rank', dimension, ranks)
                check_correlations(path, ref, actual, '_eigen0')
            elif kind == 'thermal':
                for method, name in [(1, 'tpq'), (5, 'ctpq')]:
                    path = expert_empty_rank_case(Path(str(prefix)+'_'+name), model, method=method,
                           options={'CalcMod': {'InitialVecType': 0},
                                    'ModPara': {'Lanczos_max': 5, 'initial_iv': 7, 'NumAve': 2,
                                                'LargeValue': 4 if method == 1 else 50,
                                                'ExpandCoef': 12, 'ExpecInterval': 2}})
                    executed_case_paths.append(path)
                    add_correlation_requests(path, requests, False)
                    run_case(path, probe, name, launcher, env)
                    initial, prenorms = captured_initials(path, dimension, ref.representatives, ranks)
                    expected = []
                    nsite = 3 if model == 'KondoGC' else 4
                    for sample, (state, prenorm) in enumerate(zip(initial, prenorms)):
                        check_ownership(path, f'sector_initial_sample{sample}_step0.rank', dimension, ranks)
                        history = dict(SS=[], Norm=[], Flct=[])
                        norm = prenorm
                        for step in range(NSTEPS):
                            hv = h@state
                            energy = np.vdot(state, hv).real
                            m = {key: abs(ref.basis@state)**2@ref.moments[key]
                                 for key in ('N', 'N2', 'D', 'D2', 'Sz', 'Sz2')}
                            beta = 2*step/(nsite*4-energy) if method == 1 else .02*step
                            history['SS'].append([beta, energy, np.vdot(hv, hv).real, m['D'], m['N'], step])
                            history['Norm'].append([beta, norm, prenorm, step])
                            if method == 5 or step:
                                history['Flct'].append([beta, m['N'], m['N2'], m['D'], m['D2'], m['Sz'], m['Sz2'], step])
                            if step % 2 == 0:
                                check_correlations(path, ref, state, f'_set{sample}step{step}')
                            state, norm = tpq_step(h, state, nsite) if method == 1 else ctpq_step(h, state, .02)
                        expected.append({key: np.asarray(value) for key, value in history.items()})
                    check_data(path, False, expected)
                    for sample in range(NUM_AVE):
                        for family in ('SS', 'Norm', 'Flct'):
                            actual = np.loadtxt(path/'output'/f'zvo_{family}_rand{sample}.dat', ndmin=2)
                            errors['data'] = max(errors['data'], float(np.max(abs(actual-expected[sample][family]))))
            else:
                # Start from CG, then replace its payload by a nonstationary
                # state so a missing TE application cannot pass this test.
                seed = expert_empty_rank_case(Path(str(prefix)+'_seed'), model, method=3,
                                               options={'CalcMod': {'OutputEigenVec': 1}})
                executed_case_paths.append(seed)
                run_case(seed, probe, 'seed', launcher, env)
                checkpoint_vector(seed, dimension, ranks)
                from symmetry_kondo_checkpoint import store_parts
                files = [seed/'output'/f'zvo_eigenvec_0_rank_{r}.dat' for r in range(ranks)]
                parts = [read_kondo_checkpoint(file) for file in files]
                parts = [(header, vector[int(header[14]):int(header[14]+header[15])]) for header, _ in parts]
                store_parts(files, parts, refresh=True)
                methods = (3, 4) if kind == 'checkpoint' else (4,)
                sources = [(seed, 'zvo_eigenvec_0', vector)]
                for source_index in range(2 if kind == 'checkpoint' else 1):
                    source, label, initial = sources[source_index]
                    for method in methods:
                        path = expert_empty_rank_case(Path(str(prefix)+f'_from{source_index}_to{method}'), model,
                            method=method, options={'CalcMod': {'InputEigenVec': 1, 'OutputEigenVec': 1},
                                'ModPara': {'Lanczos_max': 5, 'ExpandCoef': 8, 'ExpecInterval': 1, 'OutputInterval': 1}})
                        executed_case_paths.append(path)
                        attach_seed(path, source, label)
                        if method == 3:
                            for rank in range(ranks):
                                shutil.copy2(source/'output'/f'{label}_rank_{rank}.dat',
                                             path/'output'/f'zvo_eigenvec_0_rank_{rank}.dat')
                        add_correlation_requests(path, requests, False)
                        if method == 4:
                            schedule(path, 'TEOneBody', [], TIMES)
                        run_case(path, probe, 'import', launcher, env)
                        stem = 'sector_initial_sample0_step0.rank' if method == 4 else 'sector_final_sample0_step0.rank'
                        captured = read_vector_parts(path, stem, dimension)
                        np.testing.assert_allclose(captured, initial, atol=1e-12, rtol=0)
                        check_ownership(path, stem, dimension, ranks)
                        state, exact = initial.copy(), initial.copy()
                        history = dict(SS=[], Norm=[], Flct=[])
                        for step in range(5 if method == 4 else 1):
                            if method == 4:
                                dt = TIMES[step]-(TIMES[step-1] if step else TIMES[0])
                                state, norm = te_step(h, state, dt)
                                state /= norm
                                exact = exact_endpoint_step(h, exact, dt)
                            if method == 4:
                                hv = h@state
                                m = {key: abs(ref.basis@state)**2@ref.moments[key]
                                     for key in ('N', 'N2', 'D', 'D2', 'Sz', 'Sz2')}
                                history['SS'].append([TIMES[step], np.vdot(state, hv).real, np.vdot(hv, hv).real, m['D'], m['N'], step])
                                history['Norm'].append([TIMES[step], norm, step])
                                history['Flct'].append([TIMES[step], m['N'], m['N2'], m['D'], m['D2'], m['Sz'], m['Sz2'], step])
                            actual = checkpoint_vector(path, dimension, ranks, f'zvo_eigenvec_{step}')
                            error = np.linalg.norm(actual-state)
                            errors['vector'] = max(errors['vector'], error)
                            errors['endpoint'] = max(errors['endpoint'], np.linalg.norm(state-exact))
                            assert error < 3e-11 and np.linalg.norm(state-exact) < 1e-8
                            check_correlations(path, ref, state, f'_step{step}' if method == 4 else '_eigen0')
                            header, _ = read_kondo_checkpoint(path/'output'/f'zvo_eigenvec_{step}_rank_0.dat')
                            assert int(header[22]) == method
                        if method == 4:
                            assert np.linalg.norm(state-initial) > 1e-5, 'TE fixture must evolve nontrivially'
                            for family, wanted in history.items():
                                actual = np.loadtxt(path/'output'/f'zvo_{family}.dat', ndmin=2)
                                wanted = np.asarray(wanted)
                                assert np.isfinite(actual).all()
                                np.testing.assert_array_equal(actual[:, -1], wanted[:, -1])
                                error = float(np.max(abs(actual-wanted)))
                                assert error < 4e-11, (family, error)
                                errors['data'] = max(errors['data'], error)
                        if source_index == 0 and method == 4:
                            sources.append((path, 'zvo_eigenvec_4', state.copy()))
            check_manifests(executed_case_paths, model, layout, dimension, ranks,
                            env.get('OMP_NUM_THREADS'))
            print(f'{kind} {model} dim={dimension} np={ranks} OMP={env.get("OMP_NUM_THREADS", "default")} '
                  f'layout={layout} empty={max(0, ranks-dimension)} passed', flush=True)
    print('maximum errors:', errors, flush=True)
