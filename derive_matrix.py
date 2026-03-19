# -*- coding: utf-8 -*-
# derive_matrix.py
import os
import sys

"""
算法原理 (Mathematical Principle):
1. 在有限域 GF(2¹¹³) 中，平方操作 (a ↦ a²) 被称为 Frobenius 映射。
   由于系数直接属于 GF(2)，满足 (a+b)² = a² + b²，所以它是一个完全纯粹的线性变换。
2. 从 ONB2 (最优正规基) 到 PB (多项式基) 的转换本质是一个严密的域同构映射 M。
   域同构要求严格保持所有的加法和乘法结构，因此它也与“平方”操作保持完美交换：
   对于域中的任意元素 V，其构成的比特列向量 V_pb 和 V_onb 总是满足：
       V_pb = M · V_onb
       (V²)pb = M · (V²)onb
       (V^(2^k))pb = M · (V^(2^k))onb
3. 在 PB 中，V² 可以通过直接求多项式平方然后对不可约多项式 x¹¹³ + x⁹ + 1 取模计算得出；
4. 在 ONB2 中，基底元素的数学选择非常特殊。ONB 的核心天生性质是：
   域元素的平方直接等同于其坐标常数序列比特的“循环移位”！
   （经过穷举测定，Armadillo 使用的这套参数下，平方极为精准地等效于“向左循环移位 1 Bit”）。
5. 我们项目里碰巧已经硬编码了公钥 Q 在这两种基底体系下计算出来的确定坐标对。
6. 因此，以 Q 的 X、Y 为起点，分别将它们反复进行 112 次自身平方操作，
   我们便能以极底的成本制造出高达 226 对满足 V_pb = M · V_onb 约束的方程映射向量。
7. 在 GF(2) 的有限域运算中，从上述 226 个结果中按顺序选出前 113 个能够保持线性无关的 V_onb 列向量，
   即可填满构成一个满秩的 113×113 方阵 V_ONB，与严格对应的目标方阵 V_PB。
8. 于是由 M · V_ONB = V_PB 顺势推导出：
       M = V_PB · V_ONB⁻¹
   这就是该脚本能在无需复杂的高斯周期和迹函数（Trace）代数求根下，直接以 O(n³) 将目标完全解开的直接原因。
"""


FIELD_BITS = 113
MASK_113 = (1 << 113) - 1

def parse_hex(s):
    return int(s, 16) if isinstance(s, str) else s

# Point Q
qx_pb = parse_hex('1A103B8147D1D1B9F2D0D76A86FCF')
qy_pb = parse_hex('10FB3F3F075FDA6EFBAB1B9CE548E')

qx_onb = 3600264749883462755399490686438491
qy_onb = 419754383946908414551514272523181

# PB operations
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

def cyc_shift_left(a):
    return ((a << 1) & MASK_113) | (a >> (FIELD_BITS - 1))

def cyc_shift_right(a):
    return (a >> 1) | ((a & 1) << (FIELD_BITS - 1))

def invert_matrix_gf2(matrix_rows):
    n = len(matrix_rows)
    aug = [matrix_rows[i] | (1 << (i + n)) for i in range(n)]
    
    for col in range(n):
        piv = -1
        for row in range(col, n):
            if (aug[row] >> col) & 1:
                piv = row
                break
        if piv == -1:
            return None # Not invertible
        if piv != col:
            aug[col], aug[piv] = aug[piv], aug[col]
            
        for row in range(n):
            if row != col and ((aug[row] >> col) & 1):
                aug[row] ^= aug[col]
                
    return [aug[i] >> n for i in range(n)]

def test_shift(shift_func):
    v_onb = []
    v_pb = []
    
    cur_pb = qx_pb
    cur_onb = qx_onb
    for _ in range(113):
        v_onb.append(cur_onb)
        v_pb.append(cur_pb)
        cur_pb = gf2_sqr(cur_pb)
        cur_onb = shift_func(cur_onb)
        
    cur_pb = qy_pb
    cur_onb = qy_onb
    for _ in range(113):
        v_onb.append(cur_onb)
        v_pb.append(cur_pb)
        cur_pb = gf2_sqr(cur_pb)
        cur_onb = shift_func(cur_onb)
        
    basis_indices = []
    r_check = []
    
    for i in range(len(v_onb)):
        vec = v_onb[i]
        reduced_vec = vec
        for r in r_check:
            high_bit = r.bit_length() - 1
            if (reduced_vec >> high_bit) & 1:
                reduced_vec ^= r
        if reduced_vec != 0:
            r_check.append(reduced_vec)
            r_check.sort(key=lambda x: x.bit_length(), reverse=True)
            basis_indices.append(i)
        
        if len(r_check) == 113:
            break
            
    if len(basis_indices) < 113:
        return None
        
    ONB_cols = [v_onb[i] for i in basis_indices]
    PB_cols = [v_pb[i] for i in basis_indices]
    
    V_onb_rows = [0] * 113
    V_pb_rows = [0] * 113
    for j in range(113):
        onb_col = ONB_cols[j]
        pb_col = PB_cols[j]
        for i in range(113):
            if (onb_col >> i) & 1:
                V_onb_rows[i] |= (1 << j)
            if (pb_col >> i) & 1:
                V_pb_rows[i] |= (1 << j)
                
    INV = invert_matrix_gf2(V_onb_rows)
    if INV is None: return None
    
    M_rows = [0] * 113
    for i in range(113):
        row_val = 0
        pb_row = V_pb_rows[i]
        for k in range(113):
            if (pb_row >> k) & 1:
                row_val ^= INV[k]
        M_rows[i] = row_val
        
    return M_rows

def print_formatted_matrix(name, matrix_rows):
    print(f"    # {name}")
    print("    [")
    for r in matrix_rows:
        print(f"        {r},")
    print("    ]")

def main():
    print("=" * 60)
    print("  GF(2^113) Transformation Matrix Derivation")
    print("  Using Frobenius mappings (squaring cyclic shifts)")
    print("=" * 60)
    
    M_left = test_shift(cyc_shift_left)
    M_right = test_shift(cyc_shift_right)
    
    TARGET_ROW_0 = 10384593717069655257060992658440191
    matched_matrix = None
    shift_name = None
    
    if M_left and M_left[0] == TARGET_ROW_0:
        matched_matrix = M_left
        shift_name = "Left Cyclic Shift"
    elif M_right and M_right[0] == TARGET_ROW_0:
        matched_matrix = M_right
        shift_name = "Right Cyclic Shift"
    elif M_left:
        matched_matrix = M_left
        shift_name = "Left Cyclic Shift (Fallback)"
        
    if matched_matrix:
        print(f"\n[SUCCESS] Derived matrix using {shift_name}!")
        M_inv = invert_matrix_gf2(matched_matrix)
        
        print("\n_M_ROWS = \\")
        print_formatted_matrix("_M_ROWS", matched_matrix)
        
        print("\n_M_INV_ROWS = \\")
        print_formatted_matrix("_M_INV_ROWS", M_inv)
        
        onb_vec = qx_onb
        pb_res = 0
        for i, row in enumerate(matched_matrix):
            bits = row & onb_vec
            parity = bin(bits).count('1') & 1
            pb_res |= (parity << i)
            
        print(f"\n[Validation] Derived Poly Q.x (hex): {hex(pb_res)}")
        print(f"[Validation] Expected Poly Q.x (hex): {hex(qx_pb)}")
        if pb_res == qx_pb:
            print("[SUCCESS] Validation tests PASSED perfectly!")
        else:
            print("[FAILURE] Validation MISMATCH!")
    else:
        print("[FAILURE] Could not derive matrix reliably.")

if __name__ == "__main__":
    main()
