import math
import sys
import numpy as np

def main():
    #[s] set param.
    mag_J         = 1.0
    mag_h         = 1.0
    size          = 6

    if sys.argv[1] == "generate":
        func_generate(size, mag_J, mag_h)
    elif sys.argv[1] == "aft":
        func_aft(size, mag_J, mag_h)

def func_aft(size,mag_J,mag_h):
    #[s] ene
    file_path = 'output/zvo_energy.dat'
    with open(file_path, 'r') as file:
        for line in file:
            if 'Energy' in line:
                parts = line.split()
                ene   = float(parts[-1])  # last component
    #[e] ene
    #[s] green1.def
    file_name = "output/zvo_cisajs_eigen0.dat"
    data = np.loadtxt(file_name)
    m_x = (data[0][4]+data[1][4]+data[2][4]+data[3][4])/math.sqrt(2.0)
    m_y = (data[4][5]-data[5][5]+data[6][5]-data[7][5])/math.sqrt(2.0)
    m_z = (data[8][4]-data[9][4])
    print(m_x,m_y,m_z)
    #[e] green1.def
    #[s] green3.def
    file_name = "output/zvo_ThreeBody_eigen0.dat"
    data = np.loadtxt(file_name)
    def calculate_sum(data):
        total_1 = 0.0
        total_2 = 0.0
        for row in data:
            sign = 1
            if tuple(row[4:8]) in [(1, 0, 1, 1), (1, 1, 1, 2)]:
                sign *= -1
            if tuple(row[8:12]) == (2, 2, 2, 2):
                sign *= -1

            total_1 += sign * (row[12])
            total_2 += sign * (row[13])

        return total_1,total_2

    r1,r2  = calculate_sum(data)
    result = (1j/2.0)*(r1+1j*r2)
    #print("SxSySz:", result.real,result.imag)
    #[e] green3.def
    with open("result_hphi_size%d.dat"%(size), 'w') as file:
        file.write("      %.12f\n" %(mag_h) )
        file.write("      %.12f\n" %(ene) )
        file.write("      %.12f\n" %(m_x) )
        file.write("      %.12f\n" %(m_y) )
        file.write("      %.12f\n" %(m_z) )
        file.write("      %.12f\n" %(result.real) )
    
def func_generate(size, mag_J, mag_h):
    file_mag      = "mag.def"
    file_green1   = "green1.def"
    file_green3   = "green3.def"
    file_interall = "open_interall.def"
    file_name     = "open_namelist.def"
    #[s] green1.def
    num_green1    = 4+4+2
    with open(file_green1, 'w') as file:
        # ヘッダーの書き込み
        file.write("=====\n")
        file.write(f"NumGreen       {num_green1}\n")
        file.write("=====\n")
        file.write("=====\n")
        file.write("=====\n")
        for site in range(size):
            if site == 0: # Sx
                file.write(f"{site} 0 {site} 1 \n")
                file.write(f"{site} 1 {site} 0 \n")
                file.write(f"{site} 1 {site} 2 \n")
                file.write(f"{site} 2 {site} 1 \n")
            elif site == 1: #Sy
                file.write(f"{site} 0 {site} 1 \n")
                file.write(f"{site} 1 {site} 0 \n")
                file.write(f"{site} 1 {site} 2 \n")
                file.write(f"{site} 2 {site} 1 \n")
            elif site == 2: #Sz
                file.write(f"{site} 0 {site} 0 \n")
                file.write(f"{site} 2 {site} 2 \n")
    #[e] green1.def
 
    #[s] green3.def
    indices = {
        'i': 0,
        'j': 1,
        'k': 2
    }
    s_values = [(0, 1), (1, 2), (1, 0), (2, 1)]
    t_values = [(0, 1), (1, 2), (1, 0), (2, 1)]
    u_values = [(0, 0), (2, 2)]
    combinations = generate_combinations(s_values, t_values, u_values, indices)
    with open(file_green3, 'w') as file:
        file.write("=====\n")
        file.write(f"NumGreen       {len(combinations)}\n")
        file.write("=====\n")
        file.write("=====\n")
        file.write("=====\n")
        for comb in combinations:
            file.write(comb + '\n')
    #[e] green3.def
    #[s] mag.def
    num_mag  = 2*int(size/3)+4*2*int(size/3)
    with open(file_mag, 'w') as file:
        # ヘッダーの書き込み
        file.write("=====\n")
        file.write(f"NumMag {num_mag}     \n")
        file.write("=====\n")
        file.write("=====\n")
        file.write("=====\n")
        
        for site in range(size):
            if site%3 == 0: # Sx
                file.write(f"{site} 0 {site} 1 {-mag_h/math.sqrt(2.0)} 0.0\n")
                file.write(f"{site} 1 {site} 0 {-mag_h/math.sqrt(2.0)} 0.0\n")
                file.write(f"{site} 1 {site} 2 {-mag_h/math.sqrt(2.0)} 0.0\n")
                file.write(f"{site} 2 {site} 1 {-mag_h/math.sqrt(2.0)} 0.0\n")
            elif site%3 == 1:
                file.write(f"{site} 0 {site} 1 0.0 {mag_h/math.sqrt(2.0)}\n")
                file.write(f"{site} 1 {site} 0 0.0 {-mag_h/math.sqrt(2.0)}\n")
                file.write(f"{site} 1 {site} 2 0.0 {mag_h/math.sqrt(2.0)}\n")
                file.write(f"{site} 2 {site} 1 0.0 {-mag_h/math.sqrt(2.0)}\n")
            elif site%3 == 2:
                file.write(f"{site} 0 {site} 0 {-mag_h} 0.0\n")
                file.write(f"{site} 2 {site} 2 {mag_h} 0.0\n")
    #[e] mag.def
    #[s] open_interall.def
    num_J    = (4+4+4)*(size-1) # size-1 = open boundary condition
    with open(file_interall, 'w') as file:
        file.write("=====\n")
        file.write(f"NumJ {num_J}     \n")
        file.write("=====\n")
        file.write("=====\n")
        file.write("=====\n")
        for site in range(size-1):
            #[s] S+*S- & S-*S+
            file.write(f"{site}   0 {site}   1 {site+1} 1 {site+1} 0 {mag_J} 0.0\n")
            file.write(f"{site+1} 0 {site+1} 1 {site}   1 {site}   0 {mag_J} 0.0\n")
            file.write(f"{site}   0 {site}   1 {site+1} 2 {site+1} 1 {mag_J} 0.0\n")
            file.write(f"{site+1} 1 {site+1} 2 {site}   1 {site}   0 {mag_J} 0.0\n")
            file.write(f"{site}   1 {site}   2 {site+1} 1 {site+1} 0 {mag_J} 0.0\n")
            file.write(f"{site+1} 0 {site+1} 1 {site}   2 {site}   1 {mag_J} 0.0\n")
            file.write(f"{site}   1 {site}   2 {site+1} 2 {site+1} 1 {mag_J} 0.0\n")
            file.write(f"{site+1} 1 {site+1} 2 {site}   2 {site}   1 {mag_J} 0.0\n")
            #[e] S+*S- & S-*S+
            #[s] Sz*Sz
            file.write(f"{site} 0 {site} 0 {site+1} 0 {site+1} 0 {mag_J} 0.0\n")
            file.write(f"{site} 2 {site} 2 {site+1} 0 {site+1} 0 {-mag_J} 0.0\n")
            file.write(f"{site} 0 {site} 0 {site+1} 2 {site+1} 2 {-mag_J} 0.0\n")
            file.write(f"{site} 2 {site} 2 {site+1} 2 {site+1} 2 {mag_J} 0.0\n")
            #[e] Sz*Sz
    #[e] open_interall.def
    content = f"""\
        ModPara         modpara.def
        LocSpin         locspn.def
        Trans           {file_mag}
        InterAll        {file_interall}
        OneBodyG        {file_green1}
        TwoBodyG        greentwo.def
        ThreeBodyG      {file_green3}
        CalcMod         calcmod.def
        PairExcitation  pair.def
        SpectrumVec     zvo_eigenvec_0
"""
    with open(file_name, 'w') as file:
        file.write(content)

    stan = f"""\
        L={size}
        model = "SpinGC"
        method = "CG"
        lattice = "chain"
        J = {mag_J}
        2S  = 2
        H = {mag_h}
"""
    with open("stan.in", 'w') as file:
        file.write(stan)


def generate_combinations(s_values, t_values, u_values, indices):
    combinations = []
    for s in s_values:
        for t in t_values:
            for u in u_values:
                # 文字列フォーマットで組み合わせを生成
                combination = f"{indices['i']} {s[0]} {indices['i']} {s[1]} {indices['j']} {t[0]} {indices['j']} {t[1]} {indices['k']} {u[0]} {indices['k']} {u[1]}"
                combinations.append(combination)
    return combinations


if __name__ == "__main__":
    main()
