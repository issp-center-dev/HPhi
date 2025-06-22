# Output file name
output_file = "list_inv_temp.def"

# Fixed parameters
expand_coef   = 6
green_flag    = 0
wavefunc_flag = 0

# Create a list of temperatures
temperature_list = []

# Infinite temperature (β = 0)
temperature_list.append(float('inf'))

# High temperature region: 100, 90, ..., 10
temperature_list += list(range(100, 9, -10))

# Mid temperature region: 10, 9, ..., 2
temperature_list += list(range(10, 1, -1))

# Mid temperature (finer resolution): 1.8, 1.6, ..., 1.0
temperature_list += [round(x * 0.1, 5) for x in range(18, 9, -2)]

# Low temperature region: 0.99, 0.97, ..., 0.11
temperature_list += [round(x * 0.01, 5) for x in range(99, 9, -2)]

# Ultra-low temperature region: 0.099, 0.098, ..., 0.01
temperature_list += [round(x * 0.001, 5) for x in range(99, 9, -1)]

#print(temperature_list)

# Write inverse temperature (beta) and other parameters to the output file
with open(output_file, "w") as f:
    for T in temperature_list:
        if T == float('inf'):
            beta = 0.0
        else:
            beta = 1.0 / T
        f.write(f"{beta:.10f} {expand_coef} {green_flag} {wavefunc_flag}\n")

