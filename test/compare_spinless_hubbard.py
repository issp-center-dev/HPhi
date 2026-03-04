#!/usr/bin/env python3
"""
Compare SpinlessFermion with Hubbard (all up-spin) for Green's functions.

For Hubbard with 2Sz=N (all electrons up-spin), the system is equivalent to
SpinlessFermion. This script verifies one-body and two-body Green's function
implementations including off-diagonal terms.
"""

import os
import subprocess
import sys
import tempfile


def generate_greenone_entries(nsites):
    """Generate one-body Green's function entries including off-diagonal terms."""
    entries = []
    # Diagonal <n_i>
    for i in range(nsites):
        entries.append((i, 0, i, 0))
    # Nearest-neighbor off-diagonal <c^+_i c_j>
    for i in range(nsites):
        j = (i + 1) % nsites
        entries.append((i, 0, j, 0))
        entries.append((j, 0, i, 0))
    return entries


def generate_greentwo_entries(nsites):
    """Generate two-body Green's function entries including off-diagonal terms."""
    entries = []
    # Diagonal entries: <n_i n_j>
    for i in range(nsites):
        for j in range(i, nsites):
            entries.append((i, 0, i, 0, j, 0, j, 0))
    # Representative off-diagonal entries
    for i in range(nsites):
        j = (i + 1) % nsites
        k = (i + 2) % nsites
        l = (i + 3) % nsites
        entries.append((i, 0, j, 0, k, 0, k, 0))
        entries.append((i, 0, i, 0, k, 0, l, 0))
        entries.append((i, 0, j, 0, k, 0, l, 0))
    return entries


def write_locspn():
    with open("locspn.def", "w") as f:
        f.write("================================\n")
        f.write("NlocalSpin     0\n")
        f.write("================================\n")
        f.write("========i_0LocSpn_Sr=Sr_i=======\n")
        f.write("================================\n")


def write_trans(nsites):
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


def write_coulombinter(nsites, vval):
    with open("coulombinter.def", "w") as f:
        interactions = []
        for i in range(nsites):
            j = (i + 1) % nsites
            interactions.append((i, j, vval))
        f.write("=============================================\n")
        f.write("NCoulombInter          {}\n".format(len(interactions)))
        f.write("=============================================\n")
        f.write("================CoulombInter=================\n")
        f.write("=============================================\n")
        for entry in interactions:
            f.write("   {}     {}  {:.10f}\n".format(*entry))


def write_modpara_spinless(nsites, nelec):
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


def write_modpara_hubbard(nsites, nelec):
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


def write_calcmod(calc_model):
    with open("calcmod.def", "w") as f:
        f.write("CalcType   0\n")
        f.write("CalcModel   {}\n".format(calc_model))
        f.write("ReStart   0\n")
        f.write("CalcSpec   0\n")
        f.write("CalcEigenVec   0\n")
        f.write("InitialVecType   0\n")
        f.write("InputEigenVec   0\n")
        f.write("OutputEigenVec   0\n")


def write_greenone(nsites):
    entries = generate_greenone_entries(nsites)
    with open("greenone.def", "w") as f:
        f.write("===============================\n")
        f.write("NCisAjs         {}\n".format(len(entries)))
        f.write("===============================\n")
        f.write("======== Green functions ======\n")
        f.write("===============================\n")
        for e in entries:
            f.write("{}    {}    {}    {}\n".format(*e))


def write_greentwo(nsites):
    entries = generate_greentwo_entries(nsites)
    with open("greentwo.def", "w") as f:
        f.write("=============================================\n")
        f.write("NCisAjsCktAltDC        {}\n".format(len(entries)))
        f.write("=============================================\n")
        f.write("======== Green functions for Sqsuscep ======\n")
        f.write("=============================================\n")
        for e in entries:
            f.write("{}    {}    {}    {}    {}    {}    {}    {}\n".format(*e))


def write_namelist():
    with open("namelist.def", "w") as f:
        f.write("         ModPara  modpara.def\n")
        f.write("         CalcMod  calcmod.def\n")
        f.write("         LocSpin  locspn.def\n")
        f.write("           Trans  trans.def\n")
        f.write("    CoulombInter  coulombinter.def\n")
        f.write("        OneBodyG  greenone.def\n")
        f.write("        TwoBodyG  greentwo.def\n")


def generate_spinless_inputs(workdir, nsites, nelec, vval=0.5):
    """Generate SpinlessFermion input files."""
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    write_locspn()
    write_trans(nsites)
    write_coulombinter(nsites, vval)
    write_modpara_spinless(nsites, nelec)
    write_calcmod(7)  # SpinlessFermion
    write_greenone(nsites)
    write_greentwo(nsites)
    write_namelist()


def generate_hubbard_inputs(workdir, nsites, nelec, vval=0.5):
    """Generate Hubbard input files with all up-spin electrons (2Sz=N)."""
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)
    write_locspn()
    write_trans(nsites)
    write_coulombinter(nsites, vval)
    write_modpara_hubbard(nsites, nelec)
    write_calcmod(0)  # Hubbard
    write_greenone(nsites)
    write_greentwo(nsites)
    write_namelist()


def read_greenone(filename):
    """Read one-body Green's function output."""
    data = {}
    with open(filename, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 6:
                key = tuple(int(x) for x in parts[:4])
                value = complex(float(parts[4]), float(parts[5]))
                data[key] = value
    return data


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


def compare_data(spinless_data, hubbard_data, label):
    """Compare Spinless and Hubbard datasets and return max absolute difference."""
    max_diff = 0.0
    keys = sorted(set(spinless_data.keys()) | set(hubbard_data.keys()))
    for key in keys:
        sp_val = spinless_data.get(key, 0)
        hb_val = hubbard_data.get(key, 0)
        diff = abs(sp_val - hb_val)
        if diff > max_diff:
            max_diff = diff
            print(
                f"  [{label}] Key {key}: "
                f"Spinless={sp_val:.10f}, Hubbard={hb_val:.10f}, diff={diff:.10f}"
            )
    return max_diff


def main():
    if len(sys.argv) < 2:
        print("Usage: {} <HPhi_path>".format(sys.argv[0]))
        sys.exit(1)

    hphi_path = os.path.abspath(sys.argv[1])
    nsites = 4
    nelec = 2
    vval = 0.5

    with tempfile.TemporaryDirectory() as tmpdir:
        spinless_dir = os.path.join(tmpdir, "spinless")
        hubbard_dir = os.path.join(tmpdir, "hubbard")

        # Run SpinlessFermion
        generate_spinless_inputs(spinless_dir, nsites, nelec, vval)
        subprocess.run([hphi_path, "-e", "namelist.def"], capture_output=True, check=True)

        # Run Hubbard
        generate_hubbard_inputs(hubbard_dir, nsites, nelec, vval)
        subprocess.run([hphi_path, "-e", "namelist.def"], capture_output=True, check=True)

        # Compare one-body and two-body Green's functions
        spinless_onebody = read_greenone(os.path.join(spinless_dir, "output", "zvo_cisajs.dat"))
        hubbard_onebody = read_greenone(os.path.join(hubbard_dir, "output", "zvo_cisajs.dat"))
        spinless_twobody = read_greentwo(os.path.join(spinless_dir, "output", "zvo_cisajscktalt.dat"))
        hubbard_twobody = read_greentwo(os.path.join(hubbard_dir, "output", "zvo_cisajscktalt.dat"))

        max_diff_onebody = compare_data(spinless_onebody, hubbard_onebody, "OneBodyG")
        max_diff_twobody = compare_data(spinless_twobody, hubbard_twobody, "TwoBodyG")
        max_diff = max(max_diff_onebody, max_diff_twobody)

        print(f"\nMax one-body difference: {max_diff_onebody:.10f}")
        print(f"Max two-body difference: {max_diff_twobody:.10f}")
        if max_diff < 1e-8:
            print("PASSED: SpinlessFermion matches Hubbard (2Sz=N) for one-body/two-body G")
            sys.exit(0)

        print("FAILED: Results differ")
        sys.exit(1)


if __name__ == "__main__":
    main()
