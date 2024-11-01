from itertools import combinations
from sympy import Matrix, Poly, fraction, sqrt


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


def tetrahedral_angle_check(vs, phi, N):
    vals = []
    for v1, v2 in combinations(vs, 2):
        v1 = Matrix(v1)
        v2 = Matrix(v2)
        vals.append(
            v1.dot(v2).subs(phi, (1 + sqrt(5)) / 2).subs(N, sqrt(3)).equals(-1 / 3)
        )
    return all(vals)
