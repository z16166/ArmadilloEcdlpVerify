# ArmadilloEcdlpVerify

Python verification scripts for Armadillo ECDLP, operating over the elliptic curve E(GF(2¹¹³)): y² + xy = x³ + x² + 1.

## Files

- **`verify_final.py`**
  Final cryptographic verification script.
  - Verifies that the base point G and public key Q precisely reside on the elliptic curve.
  - Performs **Subgroup Order Verification**: Checks that n × G = O and n × Q = O (where O is the Point at Infinity) to ensure point membership in the prime-order subgroup.
  - Computes the private key scalar multiplication Q = k × G to validate the solution over Polynomial Basis.

- **`basis_convert.py`**
  Conversion tool for bidirectional mapping between **Optimal Normal Basis type 2 (ONB2)** and **Polynomial Basis (PB)** representations over GF(2¹¹³).

- **`gen_basepoint.py`**
  Simulation of the Armadillo C++ base point generation logic.
  - Faithfully reproduces the custom 32-bit PRNG and `mult()` logic from `keygen_random.cpp`.
  - Reconstructs the `ECC_RandomPoint` and `ECC_Embed` logic, including the iterative X-coordinate incrementation when quadratic roots are absent.
  - Validates the base point G = 2 × P derivation for a given `BasePointInit` seed (e.g., `0xC1F7F755`).

- **`derive_matrix.py`**
  Script to algebraically derive the rigorous 113×113 bit conversion matrices (M and M⁻¹) directly from the known coordinates of Q.
  - Generates isomorphic states using the Frobenius endomorphism (squaring) acting as a completely linear transformation constraint in GF(2¹¹³).
  - Dynamically forms linearly independent base projections internally leveraging cyclic shifts (ONB2's fundamental squaring property).
  - Completes GF(2) matrix reconstruction effortlessly using precise Gaussian elimination in O(n³) operations, avoiding complex period tracing.

## Usage

You can run each of these tools independently via Python:

```bash
# To perform final proof validity tests:
python verify_final.py

# To test manual PB <-> ONB2 coordination points checks:
python basis_convert.py

# To deduce and regenerate the exact GF(2) conversion matrices analytically:
python derive_matrix.py

# To simulate C++ base point generation from a 32-bit seed:
python gen_basepoint.py
```
