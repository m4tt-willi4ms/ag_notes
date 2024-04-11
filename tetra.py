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
print(f"{uv_prime=}")

Phi = (u - sp1 * v) * (u - sp2 * v) * (u - sp3 * v) * (u - sp4 * v)

# print(Phi.subs(u, uv_prime[0]))
Phi_prime = Phi.subs([(u, U), (v, V)]).subs([(U, uv_prime[0]), (V, uv_prime[1])])
print(Phi.radsimp().expand())
print((Phi_prime / exp(2 * pi * I / 3).rewrite(cos)).radsimp().expand())

print((Phi**3).radsimp().expand())
print((Phi_prime**3).radsimp().expand())

Psi = (
    (u - conjugate(sp1) * v)
    * (u - conjugate(sp2) * v)
    * (u - conjugate(sp3) * v)
    * (u - conjugate(sp4) * v)
)

# print(Phi.subs(u, uv_prime[0]))
Psi_prime = Psi.subs([(u, U), (v, V)]).subs([(U, uv_prime[0]), (V, uv_prime[1])])
print(Psi.radsimp().expand())
print((Psi_prime / exp(-2 * pi * I / 3).rewrite(cos)).radsimp().expand())

print((Psi**3).radsimp().expand())
print((Psi_prime**3).radsimp().expand())
