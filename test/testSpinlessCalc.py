#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test tool for SpinlessFermion calculations.
Generates expert mode input files for SpinlessFermion/SpinlessFermionGC models.
"""

import os
import subprocess
import argparse


def generate_locspn(nsites, filename="locspn.def"):
    """Generate locspn.def for SpinlessFermion (all sites are itinerant)."""
    with open(filename, "w") as f:
        f.write("================================\n")
        f.write("NlocalSpin     0\n")
        f.write("================================\n")
        f.write("========i_0LocSpn_Sr=Sr_i=======\n")
        f.write("================================\n")


def generate_trans(nsites, t=1.0, filename="trans.def"):
    """Generate trans.def for 1D chain with periodic boundary."""
    transfers = []
    for i in range(nsites):
        j = (i + 1) % nsites
        # c_i^+ c_j (spin 0 only for spinless)
        transfers.append((i, 0, j, 0, t, 0.0))
        transfers.append((j, 0, i, 0, t, 0.0))

    with open(filename, "w") as f:
        f.write("========================\n")
        f.write("NTransfer      {}\n".format(len(transfers)))
        f.write("========================\n")
        f.write("========i_j_s_tijs======\n")
        f.write("========================\n")
        for t_entry in transfers:
            f.write("{} {} {} {} {:.10f} {:.10f}\n".format(*t_entry))


def generate_coulombintra(nsites, filename="coulombintra.def"):
    """Generate empty coulombintra.def (no on-site interaction for spinless)."""
    with open(filename, "w") as f:
        f.write("=============================================\n")
        f.write("NCoulombIntra          0\n")
        f.write("=============================================\n")
        f.write("================CoulombIntra=================\n")
        f.write("=============================================\n")


def generate_modpara(nsites, nelec, model, method="Lanczos", filename="modpara.def"):
    """Generate modpara.def."""
    # LOBCG calculates 5 eigenstates by default
    exct = 5 if method == "LOBCG" else 1

    with open(filename, "w") as f:
        f.write("--------------------\n")
        f.write("Model_Parameters   0\n")
        f.write("--------------------\n")
        f.write("HPhi_Cal_Parameters\n")
        f.write("--------------------\n")
        f.write("CDataFileHead  zvo\n")
        f.write("CParaFileHead  zqp\n")
        f.write("--------------------\n")
        f.write("Nsite          {}\n".format(nsites))
        if model == "SpinlessFermion":
            f.write("Ncond          {}\n".format(nelec))
        f.write("Lanczos_max    2000\n")
        f.write("initial_iv     -1\n")
        f.write("exct           {}\n".format(exct))
        f.write("LanczosEps     14\n")
        f.write("LanczosTarget  2\n")
        f.write("LargeValue     12.0\n")
        f.write("NumAve         5\n")
        f.write("ExpecInterval  20\n")


def generate_calcmod(model, method="Lanczos", filename="calcmod.def"):
    """Generate calcmod.def."""
    # CalcType: 0=Lanczos, 3=CG (LOBCG)
    calc_type = 3 if method == "LOBCG" else 0

    with open(filename, "w") as f:
        f.write("#CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG\n")
        f.write("#CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC\n")
        f.write("#           6:HubbardNConserved, 7:SpinlessFermion, 8:SpinlessFermionGC\n")
        f.write("#ResrtVec = 0:not restart, 1:restart, 2:input first vector\n")
        f.write("#CalcSpec = 0:not calculate, 1:normal, 2:shifted Krylov\n")
        f.write("CalcType   {}\n".format(calc_type))
        # SpinlessFermion=7, SpinlessFermionGC=8 (see DefCommon.h)
        if model == "SpinlessFermion":
            f.write("CalcModel   7\n")
        else:
            f.write("CalcModel   8\n")
        f.write("ReStart   0\n")
        f.write("CalcSpec   0\n")
        f.write("CalcEigenVec   0\n")
        f.write("InitialVecType   0\n")
        f.write("InputEigenVec   0\n")
        f.write("OutputEigenVec   0\n")


def generate_greenone(nsites, filename="greenone.def"):
    """Generate greenone.def (one-body Green's function)."""
    entries = []
    for i in range(nsites):
        entries.append((i, 0, i, 0))

    with open(filename, "w") as f:
        f.write("===============================\n")
        f.write("NCisAjs         {}\n".format(len(entries)))
        f.write("===============================\n")
        f.write("======== Green functions ======\n")
        f.write("===============================\n")
        for e in entries:
            f.write("{}    {}    {}    {}\n".format(*e))


def generate_greentwo(nsites, filename="greentwo.def"):
    """Generate greentwo.def (two-body Green's function)."""
    entries = []
    for i in range(nsites):
        for j in range(i, nsites):
            entries.append((i, 0, i, 0, j, 0, j, 0))

    with open(filename, "w") as f:
        f.write("=============================================\n")
        f.write("NCisAjsCktAltDC        {}\n".format(len(entries)))
        f.write("=============================================\n")
        f.write("======== Green functions for Sqsuscep ======\n")
        f.write("=============================================\n")
        for e in entries:
            f.write("{}    {}    {}    {}    {}    {}    {}    {}\n".format(*e))


def generate_namelist(filename="namelist.def"):
    """Generate namelist.def for SpinlessFermion (no CoulombIntra needed)."""
    with open(filename, "w") as f:
        f.write("         ModPara  modpara.def\n")
        f.write("         CalcMod  calcmod.def\n")
        f.write("         LocSpin  locspn.def\n")
        f.write("           Trans  trans.def\n")
        f.write("        OneBodyG  greenone.def\n")
        f.write("        TwoBodyG  greentwo.def\n")


def main():
    parser = argparse.ArgumentParser(
        description="Test tool for SpinlessFermion calculations"
    )
    parser.add_argument('-p', '--path', default='./HPhi',
                        help='Path to HPhi executable')
    parser.add_argument('-m', '--model', default='SpinlessFermion',
                        choices=['SpinlessFermion', 'SpinlessFermionGC'],
                        help='Model type')
    parser.add_argument('-s', '--sites', type=int, default=8,
                        help='Number of sites')
    parser.add_argument('-n', '--nelec', type=int, default=None,
                        help='Number of electrons (default: sites/2)')
    parser.add_argument('-t', '--method', default='Lanczos',
                        choices=['Lanczos', 'LOBCG'],
                        help='Calculation method')
    parser.add_argument('--hopping', type=float, default=1.0,
                        help='Hopping parameter')
    parser.add_argument('-mpi', '--mpi', default='',
                        help='MPI command')

    args = parser.parse_args()

    nsites = args.sites
    nelec = args.nelec if args.nelec else nsites // 2
    model = args.model
    method = args.method

    # Generate input files
    generate_locspn(nsites)
    generate_trans(nsites, args.hopping)
    generate_modpara(nsites, nelec, model, method)
    generate_calcmod(model, method)
    generate_greenone(nsites)
    generate_greentwo(nsites)
    generate_namelist()

    # Run HPhi
    cmd = "{} {} -e namelist.def".format(args.mpi, args.path).strip()
    print("Running: {}".format(cmd))
    subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    main()
