from sympy import diff, symbols, exp, pi, sqrt, simplify

# Define the symbols
x, x0, t, V0_star = symbols("x x0  t V0")

# Define the function V(x, t)
V = (V0_star / sqrt(4 * pi * t)) * exp(-((x - x0)**2) / (4 * t) - t)

# Second derivative with respect to position (x)
dV_dx = simplify(diff(V, x, 2))

# Derivative with respect to time (t)
dV_dt = simplify(diff(V, t))

# Simplify the difference between the derivatives
difference = simplify(dV_dx - dV_dt)

# Print the results
print("V(x, t):", V)
print("Second derivative with respect to x:", dV_dx)
print("Derivative with respect to t:", dV_dt)
print("Simplified Difference:", difference)