# ==============================================================================
# SCRIPT 1: TT-COMPONENT CONSTRAINT SOLVER
# ==============================================================================
# PURPOSE: 
# To DERIVE (not verify) the relationship between the theory's parameters.
# We input a generic metric ansatz with unknown parameters (epsilon, zeta)
# and a generic action with unknown couplings (alpha, beta).
# The script finds what relationship MUST exist for the field equation to hold.
# ==============================================================================

import sympy
import time

# Initialize sympy's printing
sympy.init_printing(use_unicode=True)

print("="*70)
print("ALT FIELD EQUATION SOLVER: TT-COMPONENT")
print("="*70)

# --- 1. Define all necessary symbolic variables ---
r = sympy.Symbol('r', positive=True)
G, M, c = sympy.symbols('G M c')
# These are the UNKNOWN parameters we want to find relationships for:
epsilon, zeta, alpha, beta = sympy.symbols('epsilon zeta alpha beta')

# --- 2. Print the Inputs (The "No Cheating" Check) ---
print("\n[INPUT CONDITIONS]")
print("1. Metric Ansatz (Perturbed Schwarzschild):")
print("   A(r) = 1 - 2(GM/c^2 r) + epsilon * (GM/c^2 r)^2")
print("   B(r) = (1 - 2(GM/c^2 r) + zeta * (GM/c^2 r)^2)^-1")
print("\n2. Unknown Parameters (Unconstrained Symbols):")
print("   Geometric parameters: epsilon, zeta")
print("   Coupling parameters:  alpha, beta")
print("\n3. Goal:")
print("   Find the algebraic sum that makes the Field Equation = 0.")
print("-" * 70)

start_time = time.time()

print("\n[CALCULATION START]")

# --- 3. Define Metric & Calculate Derivatives ---
print("Step 3: Defining metric functions and calculating derivatives...")

# Define fields
rho_I = (G * M) / (c**2 * r)
A_r = 1 - 2*rho_I + epsilon*rho_I**2
B_r_inv = 1 - 2*rho_I + zeta*rho_I**2
B_r = 1 / B_r_inv
g_tt = -A_r * c**2

# Calculate derivatives
A_prime = sympy.diff(A_r, r)
B_prime = sympy.diff(B_r, r)
rho_I_prime = sympy.diff(rho_I, r)
rho_I_dprime = sympy.diff(rho_I_prime, r)

# --- 4. Calculate Geometric Tensor Components ---
print("Step 4: Computing geometric tensor components (G_tt, Box_rho)...")

# Kinetic term (d_rho)^2
d_rho_I_sq = B_r_inv * rho_I_prime**2

# Einstein Tensor Component G_tt
G_tt = (A_r * c**2 / (r**2 * B_r)) * (r * B_prime - B_r**2 + B_r)

# Christoffel Symbol Gamma^r_{tt} needed for covariant derivative
Gamma_r_tt = sympy.diff(A_r, r) / (2 * B_r)

# Covariant derivative term: nabla_t(nabla_t(rho_I)) = -Gamma^r_{tt} * rho_I'
cov_deriv_tt_term = -Gamma_r_tt * rho_I_prime

# d'Alembertian (Box) of rho_I
box_rho_I = B_r_inv * (rho_I_dprime + (2/r - B_prime/(2*B_r) + A_prime/(2*A_r)) * rho_I_prime)

# --- 5. Assemble Field Equation ---
print("Step 5: Assembling LHS (Shell Tensor) and RHS (Stress-Energy)...")

# LHS is 2 * S_tt, where S_tt is the ALT Geometric Tensor component
# S_tt = (1-a*rho^2)G_tt + 2a(rho*nabla_t(nabla_t(rho))) - 2a*g_tt(rho*Box(rho) + (d_rho_I)^2)
LHS = 2 * ( (1 - alpha*rho_I**2)*G_tt 
            + 2*alpha*(rho_I * cov_deriv_tt_term) 
            - 2*alpha*g_tt*(rho_I*box_rho_I + d_rho_I_sq) )

# RHS is beta * T_tt, where T_tt is the stress-energy tensor component for rho_I
# T_tt = -1/2 * g_tt * (d_rho_I)^2
RHS = beta * (-sympy.Rational(1, 2) * g_tt * d_rho_I_sq)

# The full equation we need to solve is LHS - RHS = 0
master_equation = LHS - RHS

# --- 6. Series Expansion ---
print("Step 6: Expanding series to isolate leading order terms...")
# We expand for large r. The leading order term should be 1/r^4.
# Using .expand() first helps sympy handle the complex expression.
series_expansion = sympy.series(master_equation.expand(), r, sympy.oo, 5).removeO()

# --- 7. Isolate and Clean the Results ---
print("Step 7: Extracting and cleaning the coefficient...")

# 1. Extract the raw coefficient including G, M, c
raw_coefficient = sympy.simplify(series_expansion.coeff(r**-4))

# 2. Define the Normalization factor (physical constants)
# All terms naturally scale with (G^2 M^2) / (2 c^2)
normalization_factor = (G**2 * M**2) / (2 * c**2)

# 3. Clean up constants to leave only the dimensionless parameters
final_algebraic_expression = sympy.simplify(raw_coefficient / normalization_factor)

# Define the constraint equation object
constraint_eq = sympy.Eq(final_algebraic_expression, 0)

end_time = time.time()

# --- 8. Display the Final Result (The Discovery) ---
print("\n" + "="*70)
print(f"CALCULATION COMPLETE (Total time: {end_time - start_time:.2f} seconds)")
print("="*70)

print("\n[OUTPUT ANALYSIS]")
print("The computer has isolated the leading-order residual term.")
print("Since we did not pre-set values for alpha/beta/epsilon/zeta,")
print("the computer returned a linear combination of them.")

print("\n--- STAGE 1: RAW PHYSICAL OUTPUT ---")
print("This is the exact term found by the calculus engine, including units:")
sympy.pprint(raw_coefficient)

print("\n--- STAGE 2: NORMALIZATION ---")
print("Since Gravity (G), Mass (M), and Speed of Light (c) are non-zero,")
print("we divide the raw output by the common physical factor:")
sympy.pprint(normalization_factor)

print("\n--- STAGE 3: DIMENSIONLESS RELATIONSHIP ---")
print("This reveals the pure algebraic constraint between the parameters:")
sympy.pprint(final_algebraic_expression)

print("\n" + "-"*40)
print("  FINAL CONSTRAINT EQUATION 1 (from TT component)")
print("-" * 40)
sympy.pprint(constraint_eq)
print("-" * 40)

print("\nCONCLUSION:")
print("For the theory to be consistent, this equation must equal zero.")
print("We have not solved the parameters yet; this is equation 1 of 3.")
print("="*70)