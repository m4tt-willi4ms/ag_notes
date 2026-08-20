from sympy import det, hessian, simplify, symbols
from utils import e_simp, reduce_multiply

e, u, v = symbols("e u v")
n = e**2 + e**3
m = e + e**4

iso_verts = []
for i in range(5):
    iso_verts.append(n * e**i)
    iso_verts.append(m * e**i)

f = u * v * reduce_multiply([u - x * v for x in iso_verts])
f = e_simp(f, e)
H = simplify(det(hessian(f, (u, v))) / (-121))
tI = (u**2 + v**2) * (u**2 - 2 * n * u * v - v**2) * (u**2 - 2 * m * u * v - v**2)
T = reduce_multiply(
    [e ** (15 - 3 * k) * tI.subs(v, v * e**k).expand() for k in range(5)]
)
T = e_simp(T, e)

I = -(H**3) / 12**3 / f**5
