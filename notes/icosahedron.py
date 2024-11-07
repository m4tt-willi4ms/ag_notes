from functools import reduce
from utils import e_simp
from sympy import I, Poly, cos, oo, simplify, sin, symbols, tan


def reduce_multiply(any_list):
    return reduce(lambda x, y: x * y, any_list)


zs = [0, oo, 1, -1, I, -I]


def rotate(z, alpha, prephase=1):
    if z is oo:
        return -cos(alpha / 2) / sin(alpha / 2)
    return (cos(alpha / 2) * prephase * z + sin(alpha / 2)) / (
        -sin(alpha / 2) * prephase * z + cos(alpha / 2)
    )


alpha, p, w, u, v, e = symbols("alpha p w u v e")


zps = [rotate(z, alpha) for z in zs]


tp = reduce_multiply([u - p * z * v for z in zps])
tp = simplify(tp.expand())
tp = tp.subs(tan(2 * alpha), w)

# Check of the above:
factors = [tp.subs(p, e**k) for k in range(5)]
result = reduce_multiply(factors)

result = Poly(result, u)


new_coeffs = []


for coeff in result.all_coeffs():
    new_coeffs.append(e_simp(coeff, e))


result = Poly.from_list(new_coeffs, u).as_expr()
result

T = result.subs(w, 2)
