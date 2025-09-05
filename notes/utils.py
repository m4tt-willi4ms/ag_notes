from functools import reduce
from itertools import combinations
from math import isclose

from sympy import I, Matrix, N, Poly, conjugate, diff, fraction, im, re, simplify, sqrt


def reduce_multiply(any_list):
    return reduce(lambda x, y: x * y, any_list)


def phi_simp(expr, phi):
    """
    Reduces expressions assuming phi**2 = phi + 1 (as is the case for the golden
    ratio)
    """
    e1, e2 = fraction(expr.factor())
    e1 = Poly(e1, phi)
    while len(e1.all_coeffs()) > 2:
        orig_coeffs = e1.all_coeffs()
        new_coeffs = e1.all_coeffs()[1:]
        new_coeffs[0] += orig_coeffs[0]
        new_coeffs[1] += orig_coeffs[0]
        e1 = Poly.from_list(new_coeffs, phi)
    e2 = Poly(e2, phi)
    while len(e2.all_coeffs()) > 2:
        orig_coeffs = e2.all_coeffs()
        new_coeffs = e2.all_coeffs()[1:]
        new_coeffs[0] += orig_coeffs[0]
        new_coeffs[1] += orig_coeffs[0]
        e2 = Poly.from_list(new_coeffs, phi)
    final_val = e1.as_expr() / e2.as_expr()
    return final_val.factor()


def n_simp(expr, N):
    """
    Reduces expressions assuming N**2 = 3 (as is the case if N = sqrt(3))
    """
    num, denom = fraction(expr)
    num = Poly(num, N)
    denom = Poly(denom.expand(), N)
    new_coeffs = [0, 0]
    for i, coeff in enumerate(reversed(num.all_coeffs())):
        if i % 2 == 0:
            new_coeffs[0] += 3 ** (i // 2) * coeff
        else:
            new_coeffs[1] += 3 ** ((i - 1) // 2) * coeff
    num = Poly.from_list(reversed(new_coeffs), N)
    new_coeffs = [0, 0]
    for i, coeff in enumerate(reversed(denom.all_coeffs())):
        if i % 2 == 0:
            new_coeffs[0] += 3 ** (i // 2) * coeff
        else:
            new_coeffs[1] += 3 ** ((i - 1) // 2) * coeff
    denom = Poly.from_list(reversed(new_coeffs), N)
    return (num / denom).as_expr()


def e_simp(expr, e, full=True, cleanup=True):
    result = Poly(expr, e)
    e_coeffs = reversed(result.all_coeffs())
    simp_coeffs = [0] * 5
    # enforce e**5 = 1
    for i, coeff in enumerate(e_coeffs):
        simp_coeffs[i % 5] += coeff
    # enforce e**4 = - e**3 - e**2 - e - 1
    if full:
        simp_coeffs = [simp_coeffs[i] - simp_coeffs[4] for i in range(4)]
    if cleanup:
        return simplify(Poly.from_list(reversed(simp_coeffs), e).as_expr())
    return Poly.from_list(reversed(simp_coeffs), e).as_expr()


def tetrahedral_angle_check(vs, phi, N):
    vals = []
    for v1, v2 in combinations(vs, 2):
        v1 = Matrix(v1)
        v2 = Matrix(v2)
        vals.append(
            v1.dot(v2).subs(phi, (1 + sqrt(5)) / 2).subs(N, sqrt(3)).equals(-1 / 3)
        )
    return all(vals)


def dot_product_from_complex_vals(z1, z2):
    r1 = 1 + z1 * conjugate(z1)
    x1 = (z1 + conjugate(z1)) / r1
    y1 = I * (z1 - conjugate(z1)) / r1
    z1 = (z1 * conjugate(z1) - 1) / r1
    r2 = 1 + z2 * conjugate(z2)
    x2 = (z2 + conjugate(z2)) / r2
    y2 = I * (z2 - conjugate(z2)) / r2
    z2 = (z2 * conjugate(z2) - 1) / r2
    return x1 * x2 + y1 * y2 + z1 * z2


def complex_tetrahedral_angle_check(z1, z2):
    dp = dot_product_from_complex_vals(z1, z2)
    dp_evalf = N(dp, 5)
    try:
        if not isclose(im(dp_evalf), 0.0, abs_tol=1e-5):
            return False
        else:
            return isclose(re(dp_evalf), -1 / 3, abs_tol=1e-4)
    except TypeError as e:
        print(f"Failed: {z1}, {z2}")
        return False

def combo_complex_tetrahedral_angle_check(zs):
    checks = []
    for z1, z2 in combinations(zs, 2):
        checks.append(complex_tetrahedral_angle_check(z1, z2))
    return all(checks)

def z_simp(expr, z_n):
    poly = Poly(expr, z_n)
    coeffs = poly.all_coeffs()
    while len(coeffs) > 2:
        top_factor = coeffs[0]
        coeffs[1] -= top_factor
        coeffs[2] += top_factor
        coeffs.pop(0)
    return Poly.from_list(coeffs, z_n).as_expr()


def k_simp(expr, k, z_n):
    poly = Poly(expr, k)
    coeffs = poly.all_coeffs()
    while len(coeffs) > 2:
        top_factor = coeffs[0]
        coeffs[2] += 3 * (3 + z_n) * top_factor
        coeffs.pop(0)
    return Poly.from_list(coeffs, k).as_expr()

def j_simp(expr, j, z_n):
    poly = Poly(expr, j)
    coeffs = poly.all_coeffs()
    while len(coeffs) > 2:
        top_factor = coeffs[0]
        coeffs[2] += (2 - z_n) * top_factor
        coeffs.pop(0)
    return Poly.from_list(coeffs, j).as_expr()


def two_simp(expr, i, full=True, cleanup=True):
    result = Poly(expr, i)
    e_coeffs = reversed(result.all_coeffs())
    simp_coeffs = [0] * 3
    # enforce i**3 = 1
    for k, coeff in enumerate(e_coeffs):
        simp_coeffs[k % 3] += coeff
    # enforce i**2 = - i - 1
    if full:
        simp_coeffs = [simp_coeffs[i] - simp_coeffs[2] for i in range(2)]
    if cleanup:
        return simplify(Poly.from_list(reversed(simp_coeffs), i).as_expr())
    return Poly.from_list(reversed(simp_coeffs), i).as_expr()

def hessian(f, u, v):
    du2 = diff(diff(f, u), u)
    dudv = diff(diff(f, u), v)
    dv2 = diff(diff(f, v), v)
    return du2 * dv2 - dudv**2

def jacobian(f, g, u, v):
    return diff(f, u) * diff(g, v) - diff(f, v) * diff(g, u)