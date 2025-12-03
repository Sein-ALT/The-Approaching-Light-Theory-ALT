# ==============================================================================
# SCRIPT 3: RR-COMPONENT CONSTRAINT SOLVER
# ==============================================================================
# PURPOSE: 
# To DERIVE the third and final algebraic constraint on the theory's parameters.
# We solve the radial-radial (rr) component of the geometric field equation
# in the weak-field limit.
# ==============================================================================

import sympy
import time

# Initialize sympy's printing
sympy.init_printing(use_unicode=True)

print("="*70)
print("ALT FIELD EQUATION SOLVER: RR-COMPONENT")
print("="*70)

# --- 1. Define all necessary symbolic variables ---
r = sympy.Symbol('r', positive=True)
G, M, c = sympy.symbols('G M c')
# These are the UNKNOWN parameters we want to find relationships for:
epsilon, zeta, alpha, beta = sympy.symbols('epsilon zeta alpha beta')

# --- 2. Print the Inputs ---
print("\n[INPUT CONDITIONS]")
print("1. Field Equation Component:")
print("   2 * S_rr = beta * T_rr")
print("\n2. Metric Ansatz (Perturbed Schwarzschild):")
print("   A(r) = 1 - 2(GM/c^2 r) + epsilon * (GM/c^2 r)^2")
print("   B(r) = (1 - 2(GM/c^2 r) + zeta * (GM/c^2 r)^2)^-1")
print("\n3. Goal:")
print("   Isolate the leading-order residual to close the system of equations.")
print("-" * 70)

start_time = time.time()

# --- 3. Define the math logic (The "Black Box") ---
print("\n[CALCULATION START]")
print("Computing tensors and derivatives...")

# Define fields
rho_I = (G * M) / (c**2 * r)
A_r = 1 - 2*rho_I + epsilon*rho_I**2
B_r_inv = 1 - 2*rho_I + zeta*rho_I**2
B_r = 1 / B_r_inv
g_rr = B_r

# Calculate derivatives
A_prime = sympy.diff(A_r, r)
B_prime = sympy.diff(B_r, r)
rho_I_prime = sympy.diff(rho_I, r)
rho_I_dprime = sympy.diff(rho_I_prime, r)

# Calculate Geometric Tensor Components
# Einstein Tensor G_rr
# Formula: G_rr = (B*A'/r - (B-1)/r^2) / A
G_rr = (B_r * A_prime / r - (B_r - 1) / r**2) / A_r

# Kinetic terms
d_rho_I_sq = B_r_inv * rho_I_prime**2

# d'Alembertian (Box) of rho_I
box_rho_I = B_r_inv * (rho_I_dprime + (2/r - B_prime/(2*B_r) + A_prime/(2*A_r)) * rho_I_prime)

# Covariant derivative nabla_r(nabla_r)
Gamma_r_rr = B_prime / (2*B_r)
cov_deriv_rr_term = rho_I_dprime - Gamma_r_rr * rho_I_prime

print("Geometric components calculated. Assembling master equation...")

# Assemble LHS (Geometric Shell Tensor S_rr)
# S_rr = (1-a*rho^2)G_rr + 2a(rho*nabla_r(nabla_r) + (nabla_r)^2) - 2a*g_rr(rho*Box + (d_rho)^2)
LHS = 2 * ( (1 - alpha*rho_I**2)*G_rr 
            + 2*alpha*(rho_I * cov_deriv_rr_term + rho_I_prime**2) 
            - 2*alpha*g_rr*(rho_I*box_rho_I + d_rho_I_sq) )

# Assemble RHS (Stress-Energy Tensor T_rr)
# T_rr = 1/2 * (rho')^2 (in static spherical symmetry)
RHS = beta * (sympy.Rational(1, 2) * rho_I_prime**2)

# The Field Equation Error
master_equation = LHS - RHS

# Expand series to find the leading order residual
print("Expanding series...")
# We expand for large r. 
series_expansion = sympy.series(master_equation.expand(), r, sympy.oo, 5).removeO()

# --- 7. Isolate and Clean the Results ---
print("Extracting leading coefficient...")
# Extract the raw coefficient of 1/r^4
raw_coefficient = sympy.simplify(series_expansion.coeff(r**-4))

# Define the Normalization factor
# In the rr component, terms scale as (GM/c^2 r^2)^2 ~ G^2 M^2 / c^4
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
print("The computer isolated the leading-order residual (1/r^4) from the RR equation.")

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
print("  FINAL CONSTRAINT EQUATION 3")
print("-" * 40)
sympy.pprint(constraint_eq)
print("-" * 40)

print("\nCONCLUSION:")
print("This is the third and final equation.")
print("We now have a system of 3 linear equations for 3 unknowns (epsilon, zeta, beta)")
print("in terms of alpha.")
print("="*70)