from itertools import combinations

from functools import reduce

from sympy import I, cot, exp, pi, sqrt, symbols

from utils import complex_tetrahedral_angle_check

e, k = symbols("e k")
w1 = I * k * (e**3 - e**2)
n = e**2 + e**3
m = e + e**4
w1h = -1 / w1
w6 = (w1 * n + 1) / (w1 - n)
w6h = (w1 * m + 1) / (w1 - m)


def reduce_multiply(any_list):
    return reduce(lambda x, y: x * y, any_list)


u, v = symbols("u v")

terms = []
roots = []
for i in range(5):
    roots.extend([w1 * e**i, w1h * e**i, w6 * e**i, w6h * e**i])
    terms.append(u - w1 * e**i * v)
    terms.append(u - w1h * e**i * v)
    terms.append(u - w6 * e**i * v)
    terms.append(u - w6h * e**i * v)

expr = reduce_multiply(terms[:4])

tetras = []
for i1, i2, i3, i4 in combinations(roots, 4):
    val1 = (
        i1.subs(e, exp(2 * I * pi / 5))
        .subs(k, 2 / sqrt(5) / (sqrt(3) - cot(pi / 5)))
        .evalf()
    )
    val2 = (
        i2.subs(e, exp(2 * I * pi / 5))
        .subs(k, 2 / sqrt(5) / (sqrt(3) - cot(pi / 5)))
        .evalf()
    )
    val3 = (
        i3.subs(e, exp(2 * I * pi / 5))
        .subs(k, 2 / sqrt(5) / (sqrt(3) - cot(pi / 5)))
        .evalf()
    )
    val4 = (
        i4.subs(e, exp(2 * I * pi / 5))
        .subs(k, 2 / sqrt(5) / (sqrt(3) - cot(pi / 5)))
        .evalf()
    )
    if all(
        map(
            lambda x: complex_tetrahedral_angle_check(x[0], x[1]),
            combinations([val1, val2, val3, val4], 2),
        )
    ):
        tetras.append([i1, i2, i3, i4])
        print([i1, i2, i3, i4])
tetras
