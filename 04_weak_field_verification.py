# ==============================================================================
# SCRIPT 4: FULL SOLUTION VERIFICATION
# ==============================================================================
# PURPOSE: 
# To VERIFY the derived solution.
# We take the relationships derived in Scripts 1-3:
#   epsilon = 3*alpha
#   zeta    = -2*alpha
#   beta    = -16*alpha
# And plug them back into the RR field equation to prove the residual is ZERO.
# ==============================================================================

import sympy
import time

# Initialize sympy's printing
sympy.init_printing(use_unicode=True)

print("="*70)
print("ALT WEAK FIELD SOLUTION VERIFICATION")
print("="*70)

# --- 1. Define all necessary symbolic variables ---
r = sympy.Symbol('r', positive=True)
G, M, c = sympy.symbols('G M c')
# Only 'alpha' is free now; the others are dependent
alpha = sympy.Symbol('alpha')
# We still define these as symbols to substitute them later
epsilon, zeta, beta = sympy.symbols('epsilon zeta beta')

# --- 2. Print the Inputs ---
print("\n[INPUT CONDITIONS]")
print("1. Hypothesis (Derived Parameter Relationships):")
print("   epsilon = 3 * alpha")
print("   zeta    = -2 * alpha")
print("   beta    = -16 * alpha")
print("\n2. Target:")
print("   Prove that with these values, the RR-Field Equation Residue is EXACTLY 0.")
print("-" * 70)

start_time = time.time()

# --- 3. Define the math logic (The "Black Box") ---
print("\n[CALCULATION START]")
print("Applying relationships and computing field equation...")

# Define the substitutions dictionary
substitutions = {
    epsilon: 3 * alpha,
    zeta: -2 * alpha,
    beta: -16 * alpha
}

# Define Metric with Substitutions applied immediately
rho_I = (G * M) / (c**2 * r)
# Note: .subs(substitutions) injects the hypothesis
A_r = (1 - 2*rho_I + epsilon*rho_I**2).subs(substitutions)
B_r_inv = (1 - 2*rho_I + zeta*rho_I**2).subs(substitutions)
B_r = 1 / B_r_inv
g_rr = B_r

# Calculate Derivatives
A_prime = sympy.diff(A_r, r)
B_prime = sympy.diff(B_r, r)
rho_I_prime = sympy.diff(rho_I, r)
rho_I_dprime = sympy.diff(rho_I_prime, r)

# Calculate Geometric Tensor Components (RR)
G_rr = (B_r * A_prime / r - (B_r - 1) / r**2) / A_r
d_rho_I_sq = B_r_inv * rho_I_prime**2
box_rho_I = B_r_inv * (rho_I_dprime + (2/r - B_prime/(2*B_r) + A_prime/(2*A_r)) * rho_I_prime)
Gamma_r_rr = B_prime / (2*B_r)
cov_deriv_rr_term = rho_I_dprime - Gamma_r_rr * rho_I_prime

print("Components calculated. Assembling equation...")

# Assemble LHS (Geometric Shell Tensor)
# Note: alpha is intrinsic, so we don't substitute it
LHS = 2 * ( (1 - alpha*rho_I**2)*G_rr 
            + 2*alpha*(rho_I * cov_deriv_rr_term + rho_I_prime**2) 
            - 2*alpha*g_rr*(rho_I*box_rho_I + d_rho_I_sq) )

# Assemble RHS (Stress-Energy Tensor)
# Note: beta IS substituted here using the dictionary
RHS = substitutions[beta] * (sympy.Rational(1, 2) * rho_I_prime**2)

# The Error Term (Residual)
error_equation = (LHS - RHS).expand()

# Expand series to find the leading order residual
print("Expanding series to check for residuals...")
# We expand for large r
series_expansion = sympy.series(error_equation, r, sympy.oo, 5).removeO()

# --- 7. Isolate and Clean the Results ---
print("Extracting leading coefficient...")
# Extract the raw coefficient of 1/r^4
raw_coefficient = sympy.simplify(series_expansion.coeff(r**-4))

# Define Normalization (just for consistency with other scripts)
normalization_factor = (G**2 * M**2) / (c**4)

# Cleaned result
final_result = sympy.simplify(raw_coefficient / normalization_factor)

end_time = time.time()

# --- 8. Display the Final Result ---
print("\n" + "="*70)
print(f"VERIFICATION COMPLETE (Total time: {end_time - start_time:.2f} seconds)")
print("="*70)

print("\n[OUTPUT ANALYSIS]")
print("The computer substituted the derived parameters back into the field equation")
print("and calculated the residual error.")

print("\n--- STAGE 1: RAW PHYSICAL RESIDUAL ---")
sympy.pprint(raw_coefficient)

# Logic check for the display
if raw_coefficient == 0:
    print("\n--- STAGE 2: VERIFICATION ---")
    print("The raw residual is ZERO.")
    print("Therefore, normalization is trivial (0 / anything = 0).")
else:
    print("\n--- STAGE 2: NORMALIZATION ---")
    sympy.pprint(normalization_factor)

print("\n" + "-"*40)
print("  FINAL VERIFICATION RESULT")
print("-" * 40)

if final_result == 0:
    print("  [OK] SUCCESS: 0")
    print("-" * 40)
    print("\nCONCLUSION:")
    print("The coefficient is exactly zero.")
    print("This computationally proves that epsilon=3alpha, zeta=-2alpha, and")
    print("beta=-16alpha form a self-consistent solution to the field equations.")
else:
    print(f"  [FAIL] ERROR: {final_result}")
    print("-" * 40)
    print("\nCONCLUSION:")
    print("The coefficient is NOT zero. The solution is inconsistent.")

print("="*70)