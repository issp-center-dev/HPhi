import subprocess
import numpy as np
import os

class CalcSpectrum:
    def __init__(self, filename, T_list, exct, eta, header="zvo", output_dir="./output"):
        self.filename = filename
        self.T_list = T_list
        self.exct = exct
        self.header = header
        self.output_dir = output_dir

    def Make_Spectrum_Input(self):
        for idx in range(self.exct):
            with open("calcmod.def") as f:
                lines = f.readlines()
            with open("calcmod_ex.def", "w") as fex:
                for line in lines:
                   words = line.split()
                   if words[0] == "CalcSpec" or words[0] == "OutputExVec":
                    continue
                   fex.write(line)
                fex.write("CalcSpec    1\n")

            with open("namelist.def") as f:
                lines = f.readlines()
            with open("namelist_ex_{}.def".format(idx), "w") as fex:
                for line in lines:
                    words = line.split()
                    if len(words) == 0:
                        continue
                    if words[0] == "CalcMod" or words[0] == "SpectrumVec":
                        continue
                    fex.write(line)
                fex.write("CalcMod calcmod_ex.def\n")
                fex.write("SpectrumVec    {}_eigenvec_{}\n".format(self.header, idx))

    def read_spectrum(self):
        spectrum_dict={}
        frequencies =[]
        for idx in range(self.exct):
            path_to_DG = os.path.join(self.output_dir, "{}_DynamicalGreen_{}.dat".format(header,idx))
            spectrum = np.loadtxt(path_to_DG)
            spectrum_dict[idx] = spectrum[:,2] + 1J*spectrum[:,3]
            if idx == 0 :
                frequencies = spectrum[:, 0] + 1J*spectrum[:, 1]
        self.spectrums_dict = spectrum_dict
        self.frequencies = frequencies
        return spectrum_dict

    def get_energies(self):
        energy_list = []
        with open(os.path.join(output_dir, "{}_energy.dat".format(header))) as f:
            lines = f.readlines()
            for line in lines:
                words = line.split()
                if len(words) != 0 and words[0] == "Energy":
                    energy_list.append(float(words[1]))
        self.energy_list = energy_list
        self.ene_min = energy_list[0]
        self.ene_max = energy_list[len(energy_list)-1]
        for T in T_list:
            eta_ene = np.exp(-(self.ene_max-self.ene_min)/T)
            print("T = {}: exp[-beta(ene_max-ene_mix)] = {}".format(T, eta_ene))
            if eta_ene > eta:
                print("Warning: At T = {}, eta_ene is larger than eta.".format(T))
        return energy_list

    def _calc_Z(self, T):
        Z = 0
        for ene in self.energy_list:
            ene_diff = ene-self.ene_min
            Z += np.exp(-ene_diff/T)
        return Z

    def get_finite_T_spectrum(self):
        self.read_spectrum()
        finite_T_spectrum_dict ={}
        for T in self.T_list:
            Z = self._calc_Z(T)
            spectrum = np.zeros_like(self.spectrums_dict[0])
            for idx in range(self.exct):
                spectrum += np.exp(-(self.energy_list[idx]-self.ene_min)/T)*self.spectrums_dict[idx]
            spectrum /= Z
            finite_T_spectrum_dict[T]=spectrum
        self.finite_T_spectrum_dict = finite_T_spectrum_dict
        return self.frequencies, finite_T_spectrum_dict

    def print_finite_T_spectrum(self, file_name = "Dynamical_Green"):
        for key, spectrum in self.finite_T_spectrum_dict.items():
            file_name_T = self.header + "_" + file_name + "_T_{}.dat".format(key)
            with open(os.path.join(self.output_dir, file_name_T), "w") as fw:
                for idx, value in enumerate(spectrum):
                    fw.write("{} {} {} {}\n".format(self.frequencies[idx].real, self.frequencies[idx].imag, value.real, value.imag))


path_hphi = "./HPhi"
output_dir = "./output"
filename = "TEST"
header = "zvo"
T_list = [0.01, 0.02, 0.03, 0.04, 0.05]
exct = 20
eta = 1e-4

print("Check Energy")
calcspectrum=CalcSpectrum(filename, T_list, exct, eta, header)
energy_list = calcspectrum.get_energies()
calcspectrum.Make_Spectrum_Input()
print("Calculate Spectrum")
for idx in range(exct):
    print("Process: {}/{}".format(idx, exct))
    input_path = "namelist_ex_{}.def".format(idx)
    cmd = "{} -e {} > std_{}.out".format(path_hphi, input_path, idx)
    subprocess.call(cmd, shell=True)
    cmd = "mv ./output/{}_DynamicalGreen.dat ./output/{}_DynamicalGreen_{}.dat".format(header, header, idx)
    subprocess.call(cmd, shell=True)

print("Calculate Finite-T spectrum ")
frequencies, finite_spectrum_list = calcspectrum.get_finite_T_spectrum()
calcspectrum.print_finite_T_spectrum()
