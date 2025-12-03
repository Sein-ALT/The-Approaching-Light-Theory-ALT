# ==============================================================================
# SCRIPT 2: SCALAR FIELD EQUATION SOLVER
# ==============================================================================
# PURPOSE: 
# To DERIVE the second algebraic constraint on the theory's parameters.
# We solve the scalar wave equation: beta * Box(rho) = 2 * alpha * rho * R
# in the weak-field limit to find the required parameter relationship.
# ==============================================================================

import sympy
import time

# Initialize sympy's printing
sympy.init_printing(use_unicode=True)

print("="*70)
print("ALT FIELD EQUATION SOLVER: SCALAR SECTOR")
print("="*70)

# --- 1. Define all necessary symbolic variables ---
r = sympy.Symbol('r', positive=True)
G, M, c = sympy.symbols('G M c')
# The UNKNOWN parameters:
epsilon, zeta, alpha, beta = sympy.symbols('epsilon zeta alpha beta')

# --- 2. Print the Inputs ---
print("\n[INPUT CONDITIONS]")
print("1. Scalar Field Equation:")
print("   beta * Box(rho) - 2 * alpha * rho * R = 0")
print("\n2. Metric Ansatz (Perturbed Schwarzschild):")
print("   A(r) = 1 - 2(GM/c^2 r) + epsilon * (GM/c^2 r)^2")
print("   B(r) = (1 - 2(GM/c^2 r) + zeta * (GM/c^2 r)^2)^-1")
print("\n3. Goal:")
print("   Isolate the leading-order residual to constrain alpha, beta, etc.")
print("-" * 70)

start_time = time.time()

# --- 3. Define the math logic (The "Black Box") ---
print("\n[CALCULATION START]")
print("Computing derivatives and curvature components...")

# Define fields
rho_I = (G * M) / (c**2 * r)
A_r = 1 - 2*rho_I + epsilon*rho_I**2
B_r_inv = 1 - 2*rho_I + zeta*rho_I**2
B_r = 1 / B_r_inv

# Calculate derivatives
A_prime = sympy.diff(A_r, r)
A_dprime = sympy.diff(A_prime, r)
B_prime = sympy.diff(B_r, r)
rho_I_prime = sympy.diff(rho_I, r)
rho_I_dprime = sympy.diff(rho_I_prime, r)

# Calculate d'Alembertian (Box) of rho_I
# Formula: Box(phi) = (1/sqrt(-g)) * d_mu(sqrt(-g) g^mu_nu d_nu phi)
# In static spherical symmetry:
box_rho_I = B_r_inv * (rho_I_dprime + (2/r - B_prime/(2*B_r) + A_prime/(2*A_r)) * rho_I_prime)

# Calculate Ricci Scalar R
term_R1 = (1/B_r) * (-A_dprime/A_r + (A_prime * B_prime)/(2*A_r*B_r) + (A_prime**2)/(2*A_r**2))
term_R2 = (2/r**2) * (1 - B_r_inv)
term_R3 = -(2/(r*B_r)) * (B_prime/B_r - A_prime/A_r)
R_scalar = term_R1 + term_R2 + term_R3

print("Components calculated. Assembling master equation...")

# Assemble the Field Equation
# Equation: beta * Box(rho) - 2 * alpha * rho * R = 0
master_equation = (beta * box_rho_I - 2 * alpha * rho_I * R_scalar).expand()

# Expand series to find the leading order residual
print("Expanding series...")
# We expand for large r. 
series_expansion = sympy.series(master_equation, r, sympy.oo, 5).removeO()

# --- 7. Isolate and Clean the Results ---
print("Extracting leading coefficient...")
# Extract the raw coefficient of 1/r^4
raw_coefficient = sympy.simplify(series_expansion.coeff(r**-4))

# Define the Normalization factor
# Scalar field terms scale as rho^2 ~ (1/c^2)^2 = 1/c^4
normalization_factor = (G**2 * M**2) / (c**4)

# Clean up constants
final_algebraic_expression = sympy.simplify(raw_coefficient / normalization_factor)

# Define the constraint equation object
constraint_eq = sympy.Eq(final_algebraic_expression, 0)

end_time = time.time()

# --- 8. Display the Final Result ---
print("\n" + "="*70)
print(f"CALCULATION COMPLETE (Total time: {end_time - start_time:.2f} seconds)")
print("="*70)

print("\n[OUTPUT ANALYSIS]")
print("The computer isolated the leading-order residual (1/r^4) from the scalar equation.")

print("\n--- STAGE 1: RAW PHYSICAL OUTPUT ---")
print("Exact term including units:")
sympy.pprint(raw_coefficient)

print("\n--- STAGE 2: NORMALIZATION ---")
print("Dividing by physical constants:")
sympy.pprint(normalization_factor)

print("\n--- STAGE 3: DIMENSIONLESS RELATIONSHIP ---")
print("The pure algebraic constraint:")
sympy.pprint(final_algebraic_expression)

print("\n" + "-"*40)
print("  FINAL CONSTRAINT EQUATION 2")
print("-" * 40)
sympy.pprint(constraint_eq)
print("-" * 40)

print("\nCONCLUSION:")
print("This provides the second required equation.")
print("Combined with Eq 1, we now have a system to solve for epsilon and zeta.")
print("="*70)