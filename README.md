# ArmadilloEcdlpVerify

Python verification scripts for Armadillo ECDLP, operating over the elliptic curve E(GF(2¹¹³)): y² + xy = x³ + x² + 1.

## Files

- **`verify_final.py`**
  Final cryptographic verification script.
  - Verifies that the base point G and public key Q precisely reside on the elliptic curve.
  - Computes the private key scalar multiplication Q = k × G to validate the solution over Polynomial Basis.

- **`basis_convert.py`**
  Conversion tool for bidirectional mapping between **Optimal Normal Basis type 2 (ONB2)** and **Polynomial Basis (PB)** representations over GF(2¹¹³).

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
```
