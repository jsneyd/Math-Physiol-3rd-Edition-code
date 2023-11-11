from sympy import symbols, Eq, solve, simplify, pretty

# Define symbolic variables
k,pe,p_i,kp,si,ci,km,se,ce,c0 = symbols("k,pe,p_i,kp,si,ci,km,se,ce,c0")

# Define the system of equations
eq1 = k*pe - k*p_i + kp*si*ci - km*p_i
eq2 = k*p_i - k*pe + kp*se*ce - km*pe
eq3 = k*ce - k*ci + km*p_i - kp*si*ci
eq4 = k*ci - k*ce + km*pe - kp*se*ce
eq5 = p_i + pe + ci + ce - c0

# Solve the system of equations
solution = solve([eq2, eq3, eq4, eq5], [ci, ce, p_i, pe])

# Extract individual values from the solution
ci_val, ce_val, pi_val, pe_val = solution[ci], solution[ce], solution[p_i], solution[pe]

# Calculate J
J = simplify(km*pi_val - kp*si*ci_val)

# Print the results
print("J:", pretty(simplify(J, 150)))
