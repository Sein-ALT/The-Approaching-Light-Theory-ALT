# ==============================================================================
# SCRIPT 5: GENERALIZED BIANCHI IDENTITY VERIFICATION
# ==============================================================================
# PURPOSE: 
# To VERIFY the energy-momentum conservation of the theory.
# In ALT, the geometry is non-minimally coupled to the scalar field.
# Therefore, the Shell Tensor (S_mn) is NOT divergence-free.
# We prove that: nabla^mu S_mn = alpha * rho * R * nabla_nu rho
# ==============================================================================

import sympy as sp
import time

# Initialize sympy's printing
sp.init_printing(use_unicode=True)

print("="*70)
print("ALT CONSERVATION LAW VERIFICATION")
print("="*70)

# --- 1. Define Symbols ---
alpha, rho, R = sp.symbols('alpha rho R')
# Abstract tensor terms (Scalar quantities representing contracted tensors)
Rmunu_dot_nrho = sp.symbols('Rmunu_dot_nrho')   # Represents: R_{mu nu} * nabla^mu rho
grad_rho_nu = sp.symbols('grad_rho_nu')         # Represents: nabla_nu rho
G2 = sp.symbols('G2')                           # Represents: nabla^mu rho * nabla_mu nabla_nu rho
grad_Box = sp.symbols('grad_Box')               # Represents: nabla_nu (Box rho)
Box_rho = sp.symbols('Box_rho')                 # Represents: Box rho

# --- 2. Print Inputs ---
print("\n[INPUT CONDITIONS]")
print("1. Goal:")
print("   Calculate the divergence of the Shell Tensor: nabla^mu S_{mu nu}")
print("\n2. Component Groups (Derived in Appendix J):")
print("   Group 1: From prefactor G_mn(1 - alpha*rho^2)")
print("   Group 2: From second-derivative interaction terms")
print("   Group 3: From the trace term")
print("-" * 70)

start_time = time.time()

# --- 3. Define the Math Logic ---
print("\n[CALCULATION START]")
print("Constructing divergence components...")

# Group 1: Prefactor term divergence
# nabla^mu [ (1 - alpha rho^2) G_{mu nu} ]
prefactor = -2*alpha * rho * Rmunu_dot_nrho + alpha * rho * R * grad_rho_nu

# Group 2: Divergence of second-derivative terms
# nabla^mu [ 2alpha (rho nabla_mu nabla_nu rho + nabla_mu rho nabla_nu rho) ]
second_group = 2*alpha * (2*G2 + rho*grad_Box + rho*Rmunu_dot_nrho + Box_rho*grad_rho_nu)

# Group 3: Divergence of trace term
# nabla^mu [ -2alpha g_{mu nu} (rho Box rho + (d_rho)^2) ]
trace_term = -2*alpha * (grad_rho_nu*Box_rho + rho*grad_Box + 2*G2)

# --- 4. Intermediate Expansions (Transparency) ---
print("\n--- INTERMEDIATE EXPANSIONS ---")
print("Group 1 (Prefactor):")
sp.pprint(sp.expand(prefactor))

print("\nGroup 2 (Interaction):")
sp.pprint(sp.expand(second_group))

print("\nGroup 3 (Trace):")
sp.pprint(sp.expand(trace_term))

# --- 5. Assemble and Simplify ---
print("\nSumming components...")
full_divergence = prefactor + second_group + trace_term

print("Simplifying algebraic cancellations...")
calculated_divergence = sp.simplify(full_divergence)

# --- 6. Verify Against Expectation ---
# The scalar field equation implies the exchange term must be: alpha * rho * R * grad_rho
expected_exchange = alpha * rho * R * grad_rho_nu

# Calculate Residual
residual = sp.simplify(calculated_divergence - expected_exchange)

end_time = time.time()

# --- 7. Display Final Result ---
print("\n" + "="*70)
print(f"VERIFICATION COMPLETE (Total time: {end_time - start_time:.2f} seconds)")
print("="*70)

print("\n[OUTPUT ANALYSIS]")

print("\n--- STAGE 1: CALCULATED GEOMETRIC DIVERGENCE ---")
print("The net momentum lost by the geometry (nabla^mu S_mn):")
sp.pprint(calculated_divergence)

print("\n--- STAGE 2: REQUIRED SCALAR EXCHANGE ---")
print("The momentum gained by the scalar field (from field equations):")
sp.pprint(expected_exchange)

print("\n" + "-"*40)
print("  FINAL VERIFICATION RESULT")
print("-" * 40)

if residual == 0:
    print("  [OK] SUCCESS: Residual is 0")
    print("-" * 40)
    print("\nCONCLUSION:")
    print("The geometric divergence EXACTLY matches the scalar exchange term.")
    print("Energy-Momentum is strictly conserved in the coupled system.")
else:
    print("  [FAIL] ERROR: Residual is NOT 0")
    print("-" * 40)
    print("Leftover terms:")
    sp.pprint(residual)

print("="*70)