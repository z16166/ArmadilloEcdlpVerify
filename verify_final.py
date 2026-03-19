import basis_convert

# Verification Script: Final Proof of (k, Q) validity
FIELD_BITS = 113
FIELD_MOD = (1 << 113) | (1 << 9) | 1
MASK_113 = (1 << 113) - 1

def gf2_add(a, b): return a ^ b

def gf2_reduce(a):
    while a >> FIELD_BITS:
        hi = a >> FIELD_BITS
        a = (a & MASK_113) ^ hi ^ (hi << 9)
    return a

def gf2_mul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1
        a <<= 1
        if a >> FIELD_BITS: a = (a & MASK_113) ^ (a >> FIELD_BITS) ^ ((a >> FIELD_BITS) << 9)
    return gf2_reduce(r)

def gf2_sqr(a): return gf2_mul(a, a)

def gf2_inv(a):
    if a == 0: return 0
    r = a
    for _ in range(111): r = gf2_mul(gf2_sqr(r), a)
    return gf2_sqr(r)

def ec_add(p, q):
    if p is None: return q
    if q is None: return p
    x1, y1 = p; x2, y2 = q
    if x1 == x2:
        if y1 == (x2 ^ y2): return None
        if y1 == y2: return ec_dbl(p)
        return None
    lam = gf2_mul(y1 ^ y2, gf2_inv(x1 ^ x2))
    x3 = gf2_reduce(gf2_sqr(lam) ^ lam ^ x1 ^ x2 ^ 1)
    y3 = gf2_reduce(gf2_mul(lam, x1 ^ x3) ^ x3 ^ y1)
    return (x3, y3)

def ec_dbl(p):
    if p is None: return None
    x, y = p
    if x == 0: return None
    lam = x ^ gf2_mul(y, gf2_inv(x))
    x3 = gf2_reduce(gf2_sqr(lam) ^ lam ^ 1)
    y3 = gf2_reduce(gf2_sqr(x) ^ gf2_mul(lam ^ 1, x3))
    return (x3, y3)

def ec_mul(k, p):
    if k == 0 or p is None: return None
    result = None; addend = p
    while k:
        if k & 1: result = ec_add(result, addend)
        addend = ec_dbl(addend)
        k >>= 1
    return result

def is_on_curve(p):
    if p is None: return False
    x, y = p
    lhs = gf2_sqr(y) ^ gf2_mul(x, y)
    rhs = gf2_mul(gf2_sqr(x), x) ^ gf2_sqr(x) ^ 1
    return lhs == rhs

def main():
    # Parameters for Encryptionizer
    gx_pb = int('10F2C397CA2D4ABB13FBA7BFFFA95', 16)
    gy_pb = int('1A3D81BD511A62CE29B358D2AFCCE', 16)

    qx_pb = int('1A103B8147D1D1B9F2D0D76A86FCF', 16)
    qy_pb = int('10FB3F3F075FDA6EFBAB1B9CE548E', 16)

    qx_onb2 = 3600264749883462755399490686438491
    qy_onb2 = 419754383946908414551514272523181

    k = 0x2082DF1821D6E82CEE6880211228

    qx_pb_calc = basis_convert.onb2_to_poly(qx_onb2)
    qy_pb_calc = basis_convert.onb2_to_poly(qy_onb2)
    if (qx_pb_calc == qx_pb and qy_pb_calc == qy_pb):
        print("[SUCCESS] Q in ONB2 matches Q in Polynomial Basis!")
    else:
        print("[FAILURE] Q in ONB2 does not match Q in Polynomial Basis.")
        print(hex(qx_pb_calc), hex(qx_pb))
        print(hex(qy_pb_calc), hex(qy_pb))

    qx_onb2_calc = basis_convert.poly_to_onb2(qx_pb)
    qy_onb2_calc = basis_convert.poly_to_onb2(qy_pb)
    if (qx_onb2_calc == qx_onb2 and qy_onb2_calc == qy_onb2):
        print("[SUCCESS] Q in Polynomial Basis matches Q in ONB2!")
    else:
        print("[FAILURE] Q in Polynomial Basis does not match Q in ONB2.")
        print(hex(qx_onb2_calc), hex(qx_onb2))
        print(hex(qy_onb2_calc), hex(qy_onb2))

    print("Verifying if G is on the curve...")
    if is_on_curve((gx_pb, gy_pb)):
        print("[SUCCESS] G is on the elliptic curve!")
    else:
        print("[FAILURE] G is NOT on the elliptic curve.")

    print("Verifying if Q is on the curve...")
    if is_on_curve((qx_pb, qy_pb)):
        print("[SUCCESS] Q is on the elliptic curve!")
    else:
        print("[FAILURE] Q is NOT on the elliptic curve.")

    print("Verifying k * G = Q in Encryptionizer (Polynomial Basis)...")
    Q_calc = ec_mul(k, (gx_pb, gy_pb))
    if Q_calc and Q_calc[0] == qx_pb and Q_calc[1] == qy_pb:
        print("[SUCCESS] Private key k matches Public key Q perfectly!")
    else:
        print("[FAILURE] Keys do not match.")

if __name__ == "__main__":
    main()
