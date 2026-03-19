# -*- coding: utf-8 -*-
# gen_basepoint.py

import basis_convert

FIELD_BITS = 113
MASK_113 = (1 << 113) - 1

# --- PB operations ---
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
        if a >> FIELD_BITS: 
            a = (a & MASK_113) ^ (a >> FIELD_BITS) ^ ((a >> FIELD_BITS) << 9)
    return gf2_reduce(r)

def gf2_sqr(a): 
    return gf2_mul(a, a)

def gf2_inv(a):
    if a == 0: return 0
    r = a
    for _ in range(111): 
        r = gf2_mul(gf2_sqr(r), a)
    return gf2_sqr(r)

# --- ONB2 operations (via isomorphism) ---
def onb2_multiply(a_onb, b_onb):
    a_pb = basis_convert.onb2_to_poly(a_onb)
    b_pb = basis_convert.onb2_to_poly(b_onb)
    return basis_convert.poly_to_onb2(gf2_mul(a_pb, b_pb))

def onb2_inverse(a_onb):
    a_pb = basis_convert.onb2_to_poly(a_onb)
    return basis_convert.poly_to_onb2(gf2_inv(a_pb))

def onb2_rotate_left(a_onb):
    return ((a_onb << 1) & MASK_113) | (a_onb >> 112)

def onb2_rotate_right(a_onb):
    return (a_onb >> 1) | ((a_onb & 1) << 112)

# --- Custom PRNG from keygen_random.cpp ---
m_rand = 100000000
m1 = 10000
b_rand = 31415821
a_seed = 0

def to_ulong(x):
    return x & 0xFFFFFFFF

def c_long_div(a, b):
    # C99 truncates toward zero
    res = abs(a) // abs(b)
    if (a < 0) ^ (b < 0): 
        return -res
    return res

def c_long_mod(a, b):
    # C99 modulo takes sign of dividend
    res = abs(a) % abs(b)
    if a < 0: 
        return -res
    return res

def mult_rand(p, q):
    # Replicate C++: unsigned long mult(long p, long q)
    def sign_ext(x):
        x = x & 0xFFFFFFFF
        return x - 0x100000000 if x >= 0x80000000 else x
        
    p_signed = sign_ext(p)
    q_signed = sign_ext(q)
    
    # unsigned long p1 = p / m1, p0 = p % m1, q1 = q / m1, q0 = q % m1;
    p1 = to_ulong(c_long_div(p_signed, m1))
    p0 = to_ulong(c_long_mod(p_signed, m1))
    q1 = to_ulong(c_long_div(q_signed, m1))
    q0 = to_ulong(c_long_mod(q_signed, m1))
    
    # (((p0 * q1 + p1 * q0) % m1) * m1 + p0 * q0) % m;
    # All arithmetic happens as 32-bit unsigned long!
    t1 = to_ulong(p0 * q1)
    t2 = to_ulong(p1 * q0)
    t3 = to_ulong(t1 + t2)
    t4 = to_ulong(t3 % m1)
    t5 = to_ulong(t4 * m1)
    t6 = to_ulong(p0 * q0)
    t7 = to_ulong(t5 + t6)
    
    return to_ulong(t7 % m_rand)

def InitRandomGenerator(seed):
    global a_seed
    a_seed = seed

def NextRandomRange(r):
    global a_seed
    a_seed = (mult_rand(a_seed, b_rand) + 1) % m_rand
    return ((a_seed // m1) * r) // m1

def NextRandomNumber():
    n1 = NextRandomRange(256)
    n2 = NextRandomRange(256)
    n3 = NextRandomRange(256)
    n4 = NextRandomRange(256)
    return (n1 << 24) | (n2 << 16) | (n3 << 8) | n4

# --- ECC Array Conversion ---
def to_int(e):
    return (e[0] << 96) | (e[1] << 64) | (e[2] << 32) | e[3]

def to_e(v):
    return [(v >> 96) & 0x1FFFF, (v >> 64) & 0xFFFFFFFF, (v >> 32) & 0xFFFFFFFF, v & 0xFFFFFFFF]

# --- ECC_Quadratic Simulation ---
def ECC_Quadratic(a_int, b_int):
    if a_int == 0:
        y0 = onb2_rotate_right(b_int)
        return y0, y0

    a2_int = onb2_rotate_left(onb2_inverse(a_int))
    k_int = onb2_rotate_right(onb2_multiply(b_int, a2_int))

    # Check for solution
    if (bin(k_int).count('1') & 1) != 0:
        return None

    k_e = to_e(k_int)
    x_e = [0, 0, 0, 0]

    mask = 1
    for bits in range(113):
        i = 3 - (bits // 32)
        l = 3 - ((bits + 1) // 32)

        r_bit = k_e[i] & mask
        t_bit = x_e[i] & mask
        r_bit ^= t_bit

        if l == i:
            r_bit <<= 1
            x_e[l] |= r_bit
            mask <<= 1
        else:
            mask = 1
            if r_bit:
                x_e[l] |= 1

    x_int = to_int(x_e)
    y0 = onb2_multiply(a_int, x_int)
    y1 = y0 ^ a_int
    return y0, y1

# --- ECC_PointDouble in ONB2 ---
def ECC_PointDouble(p1_x, p1_y):
    x1 = onb2_inverse(p1_x)
    y1 = onb2_multiply(x1, p1_y)
    theta = p1_x ^ y1

    theta2 = onb2_rotate_left(theta)

    p3_x = theta ^ theta2 ^ MASK_113

    y1_new = MASK_113 ^ theta
    t1 = onb2_multiply(y1_new, p3_x)
    x1_new = onb2_rotate_left(p1_x)
    p3_y = x1_new ^ t1
    return p3_x, p3_y

# --- Main Base Point Generator ---
def gen_basepoint(seed):
    InitRandomGenerator(seed)
    
    # Generate 4 words for X
    x_words = [0]*4
    for i in range(4):
        x_words[i] = NextRandomNumber()
    x_words[0] &= 0x1FFFF
    
    root_bit = x_words[3] & 1
    
    while True:
        x_onb = to_int(x_words)
        
        x2_onb = onb2_rotate_left(x_onb)
        x3_onb = onb2_multiply(x_onb, x2_onb)
        
        # f(x) = x^3 + a2*x^2 + a6 -> a2 and a6 are 1 (MASK_113 in ONB2)
        f_onb = x3_onb ^ x2_onb ^ MASK_113
        
        ys = ECC_Quadratic(x_onb, f_onb)
        if ys is not None:
            temp_x = x_onb
            temp_y = ys[root_bit]
            break
            
        # Increment LSB word as a normal 32-bit integer
        x_words[3] = (x_words[3] + 1) & 0xFFFFFFFF
        
    G_x, G_y = ECC_PointDouble(temp_x, temp_y)
    
    # Return as Polynomial Basis
    return basis_convert.onb2_to_poly(G_x), basis_convert.onb2_to_poly(G_y)

def main():
    target_seed = 0xC1F7F755
    expected_gx = int('10F2C397CA2D4ABB13FBA7BFFFA95', 16)
    expected_gy = int('1A3D81BD511A62CE29B358D2AFCCE', 16)
    
    print(f"Generating Base Point G from seed: {hex(target_seed)}")
    
    gx_pb, gy_pb = gen_basepoint(target_seed)
    
    print(f"\\nGenerated G_x (PB): {hex(gx_pb).upper().replace('0X', '0x')}")
    print(f"Generated G_y (PB): {hex(gy_pb).upper().replace('0X', '0x')}")
    print(f"Expected G_x (PB):  0x{expected_gx:X}")
    print(f"Expected G_y (PB):  0x{expected_gy:X}")
    
    if gx_pb == expected_gx and gy_pb == expected_gy:
        print("\\n[SUCCESS] Base point exactly matches the provided test case!")
    else:
        print("\\n[FAILURE] Base point mismatch! Algorithm does not align.")

if __name__ == "__main__":
    main()
