from sympy import I, cos, pi, symbols, sqrt, Rational, Matrix, exp, conjugate

u, v = symbols("u v")
U, V = symbols("U V")


sp1 = (1 + I) / (sqrt(3) - 1)
sp2 = (1 - I) / (sqrt(3) + 1)
sp3 = (-1 + I) / (sqrt(3) + 1)
sp4 = -(1 + I) / (sqrt(3) - 1)

orth = Rational(1, 2) * Matrix([[1 + I, -1 + I], [1 + I, 1 - I]])
uv = Matrix([u, v])

uv_prime = orth * uv

Phi = (u - sp1 * v) * (u - sp2 * v) * (u - sp3 * v) * (u - sp4 * v)

Phi_prime = Phi.subs([(u, U), (v, V)]).subs([(U, uv_prime[0]), (V, uv_prime[1])])

Psi = (
    (u - conjugate(sp1) * v)
    * (u - conjugate(sp2) * v)
    * (u - conjugate(sp3) * v)
    * (u - conjugate(sp4) * v)
)

Psi_prime = Psi.subs([(u, U), (v, V)]).subs([(U, uv_prime[0]), (V, uv_prime[1])])

t = u * v * (u - v) * (u + v) * (u**2 + v**2)
