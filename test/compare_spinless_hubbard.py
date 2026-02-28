#!/usr/bin/env python3
"""
Compare SpinlessFermion with Hubbard (all up-spin) for off-diagonal Green's function.

For Hubbard with 2Sz=N (all electrons up-spin), the system is equivalent to
SpinlessFermion. This script verifies the off-diagonal two-body Green's function
implementation.
"""

import os
import subprocess
import sys
import tempfile
import numpy as np


def generate_spinless_inputs(workdir, nsites, nelec, V=0.5):
    """Generate SpinlessFermion input files."""
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    # locspn.def
    with open("locspn.def", "w") as f:
        f.write("================================\n")
        f.write("NlocalSpin     0\n")
        f.write("================================\n")
        f.write("========i_0LocSpn_Sr=Sr_i=======\n")
        f.write("================================\n")

    # trans.def
    with open("trans.def", "w") as f:
        transfers = []
        for i in range(nsites):
            j = (i + 1) % nsites
            transfers.append((i, 0, j, 0, 1.0, 0.0))
            transfers.append((j, 0, i, 0, 1.0, 0.0))
        f.write("========================\n")
        f.write("NTransfer      {}\n".format(len(transfers)))
        f.write("========================\n")
        f.write("========i_j_s_tijs======\n")
        f.write("========================\n")
        for t in transfers:
            f.write("{} {} {} {} {:.10f} {:.10f}\n".format(*t))

    # coulombinter.def
    with open("coulombinter.def", "w") as f:
        interactions = []
        for i in range(nsites):
            j = (i + 1) % nsites
            interactions.append((i, j, V))
        f.write("=============================================\n")
        f.write("NCoulombInter          {}\n".format(len(interactions)))
        f.write("=============================================\n")
        f.write("================CoulombInter=================\n")
        f.write("=============================================\n")
        for entry in interactions:
            f.write("   {}     {}  {:.10f}\n".format(*entry))

    # modpara.def
    with open("modpara.def", "w") as f:
        f.write("--------------------\n")
        f.write("Model_Parameters   0\n")
        f.write("--------------------\n")
        f.write("HPhi_Cal_Parameters\n")
        f.write("--------------------\n")
        f.write("CDataFileHead  zvo\n")
        f.write("CParaFileHead  zqp\n")
        f.write("--------------------\n")
        f.write("Nsite          {}\n".format(nsites))
        f.write("Ncond          {}\n".format(nelec))
        f.write("Lanczos_max    2000\n")
        f.write("initial_iv     1\n")
        f.write("exct           1\n")
        f.write("LanczosEps     14\n")
        f.write("LanczosTarget  2\n")
        f.write("LargeValue     12.0\n")
        f.write("NumAve         5\n")
        f.write("ExpecInterval  20\n")

    # calcmod.def
    with open("calcmod.def", "w") as f:
        f.write("CalcType   0\n")
        f.write("CalcModel   7\n")  # SpinlessFermion
        f.write("ReStart   0\n")
        f.write("CalcSpec   0\n")
        f.write("CalcEigenVec   0\n")
        f.write("InitialVecType   0\n")
        f.write("InputEigenVec   0\n")
        f.write("OutputEigenVec   0\n")

    # greenone.def
    with open("greenone.def", "w") as f:
        f.write("===============================\n")
        f.write("NCisAjs         {}\n".format(nsites))
        f.write("===============================\n")
        f.write("======== Green functions ======\n")
        f.write("===============================\n")
        for i in range(nsites):
            f.write("{}    0    {}    0\n".format(i, i))

    # greentwo.def with off-diagonal entries
    with open("greentwo.def", "w") as f:
        entries = []
        # Diagonal
        for i in range(nsites):
            for j in range(i, nsites):
                entries.append((i, 0, i, 0, j, 0, j, 0))
        # Off-diagonal
        for i in range(nsites):
            j = (i + 1) % nsites
            k = (i + 2) % nsites
            l = (i + 3) % nsites
            entries.append((i, 0, j, 0, k, 0, k, 0))
            entries.append((i, 0, i, 0, k, 0, l, 0))
            entries.append((i, 0, j, 0, k, 0, l, 0))
        f.write("=============================================\n")
        f.write("NCisAjsCktAltDC        {}\n".format(len(entries)))
        f.write("=============================================\n")
        f.write("======== Green functions for Sqsuscep ======\n")
        f.write("=============================================\n")
        for e in entries:
            f.write("{}    {}    {}    {}    {}    {}    {}    {}\n".format(*e))

    # namelist.def
    with open("namelist.def", "w") as f:
        f.write("         ModPara  modpara.def\n")
        f.write("         CalcMod  calcmod.def\n")
        f.write("         LocSpin  locspn.def\n")
        f.write("           Trans  trans.def\n")
        f.write("    CoulombInter  coulombinter.def\n")
        f.write("        OneBodyG  greenone.def\n")
        f.write("        TwoBodyG  greentwo.def\n")


def generate_hubbard_inputs(workdir, nsites, nelec, V=0.5):
    """Generate Hubbard input files with all up-spin electrons (2Sz=N)."""
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    # locspn.def
    with open("locspn.def", "w") as f:
        f.write("================================\n")
        f.write("NlocalSpin     0\n")
        f.write("================================\n")
        f.write("========i_0LocSpn_Sr=Sr_i=======\n")
        f.write("================================\n")

    # trans.def - only up-spin hopping
    with open("trans.def", "w") as f:
        transfers = []
        for i in range(nsites):
            j = (i + 1) % nsites
            # Up-spin only
            transfers.append((i, 0, j, 0, 1.0, 0.0))
            transfers.append((j, 0, i, 0, 1.0, 0.0))
        f.write("========================\n")
        f.write("NTransfer      {}\n".format(len(transfers)))
        f.write("========================\n")
        f.write("========i_j_s_tijs======\n")
        f.write("========================\n")
        for t in transfers:
            f.write("{} {} {} {} {:.10f} {:.10f}\n".format(*t))

    # coulombinter.def - same as spinless
    with open("coulombinter.def", "w") as f:
        interactions = []
        for i in range(nsites):
            j = (i + 1) % nsites
            interactions.append((i, j, V))
        f.write("=============================================\n")
        f.write("NCoulombInter          {}\n".format(len(interactions)))
        f.write("=============================================\n")
        f.write("================CoulombInter=================\n")
        f.write("=============================================\n")
        for entry in interactions:
            f.write("   {}     {}  {:.10f}\n".format(*entry))

    # modpara.def
    with open("modpara.def", "w") as f:
        f.write("--------------------\n")
        f.write("Model_Parameters   0\n")
        f.write("--------------------\n")
        f.write("HPhi_Cal_Parameters\n")
        f.write("--------------------\n")
        f.write("CDataFileHead  zvo\n")
        f.write("CParaFileHead  zqp\n")
        f.write("--------------------\n")
        f.write("Nsite          {}\n".format(nsites))
        f.write("Ncond          {}\n".format(nelec))
        f.write("2Sz            {}\n".format(nelec))  # All up-spin
        f.write("Lanczos_max    2000\n")
        f.write("initial_iv     1\n")
        f.write("exct           1\n")
        f.write("LanczosEps     14\n")
        f.write("LanczosTarget  2\n")
        f.write("LargeValue     12.0\n")
        f.write("NumAve         5\n")
        f.write("ExpecInterval  20\n")

    # calcmod.def
    with open("calcmod.def", "w") as f:
        f.write("CalcType   0\n")
        f.write("CalcModel   0\n")  # Hubbard
        f.write("ReStart   0\n")
        f.write("CalcSpec   0\n")
        f.write("CalcEigenVec   0\n")
        f.write("InitialVecType   0\n")
        f.write("InputEigenVec   0\n")
        f.write("OutputEigenVec   0\n")

    # greenone.def - up-spin only
    with open("greenone.def", "w") as f:
        f.write("===============================\n")
        f.write("NCisAjs         {}\n".format(nsites))
        f.write("===============================\n")
        f.write("======== Green functions ======\n")
        f.write("===============================\n")
        for i in range(nsites):
            f.write("{}    0    {}    0\n".format(i, i))

    # greentwo.def with off-diagonal entries (up-spin only)
    with open("greentwo.def", "w") as f:
        entries = []
        # Diagonal
        for i in range(nsites):
            for j in range(i, nsites):
                entries.append((i, 0, i, 0, j, 0, j, 0))
        # Off-diagonal
        for i in range(nsites):
            j = (i + 1) % nsites
            k = (i + 2) % nsites
            l = (i + 3) % nsites
            entries.append((i, 0, j, 0, k, 0, k, 0))
            entries.append((i, 0, i, 0, k, 0, l, 0))
            entries.append((i, 0, j, 0, k, 0, l, 0))
        f.write("=============================================\n")
        f.write("NCisAjsCktAltDC        {}\n".format(len(entries)))
        f.write("=============================================\n")
        f.write("======== Green functions for Sqsuscep ======\n")
        f.write("=============================================\n")
        for e in entries:
            f.write("{}    {}    {}    {}    {}    {}    {}    {}\n".format(*e))

    # namelist.def
    with open("namelist.def", "w") as f:
        f.write("         ModPara  modpara.def\n")
        f.write("         CalcMod  calcmod.def\n")
        f.write("         LocSpin  locspn.def\n")
        f.write("           Trans  trans.def\n")
        f.write("    CoulombInter  coulombinter.def\n")
        f.write("        OneBodyG  greenone.def\n")
        f.write("        TwoBodyG  greentwo.def\n")


def read_greentwo(filename):
    """Read two-body Green's function output."""
    data = {}
    with open(filename, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 10:
                key = tuple(int(x) for x in parts[:8])
                value = complex(float(parts[8]), float(parts[9]))
                data[key] = value
    return data


def main():
    if len(sys.argv) < 2:
        print("Usage: {} <HPhi_path>".format(sys.argv[0]))
        sys.exit(1)

    hphi_path = sys.argv[1]
    nsites = 4
    nelec = 2
    V = 0.5

    with tempfile.TemporaryDirectory() as tmpdir:
        spinless_dir = os.path.join(tmpdir, "spinless")
        hubbard_dir = os.path.join(tmpdir, "hubbard")

        # Run SpinlessFermion
        generate_spinless_inputs(spinless_dir, nsites, nelec, V)
        subprocess.run([hphi_path, "-e", "namelist.def"], capture_output=True)

        # Run Hubbard
        generate_hubbard_inputs(hubbard_dir, nsites, nelec, V)
        subprocess.run([hphi_path, "-e", "namelist.def"], capture_output=True)

        # Compare results
        spinless_data = read_greentwo(os.path.join(spinless_dir, "output", "zvo_cisajscktalt.dat"))
        hubbard_data = read_greentwo(os.path.join(hubbard_dir, "output", "zvo_cisajscktalt.dat"))

        max_diff = 0.0
        for key in spinless_data:
            sp_val = spinless_data[key]
            hb_val = hubbard_data.get(key, 0)
            diff = abs(sp_val - hb_val)
            if diff > max_diff:
                max_diff = diff
                print(f"  Key {key}: Spinless={sp_val:.10f}, Hubbard={hb_val:.10f}, diff={diff:.10f}")

        print(f"\nMax difference: {max_diff:.10f}")
        if max_diff < 1e-8:
            print("PASSED: SpinlessFermion matches Hubbard (2Sz=N)")
            sys.exit(0)
        else:
            print("FAILED: Results differ")
            sys.exit(1)


if __name__ == "__main__":
    main()
