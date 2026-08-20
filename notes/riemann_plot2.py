import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import hsv_to_rgb


# -----------------------------
# 1) Choose your complex function
# -----------------------------
def H(u, v):
    return u**20 - 228 * (u**15 * v**5 - u**5 * v**15) + 494 * u**10 * v**10 + v**20


def f(u, v):
    return u * v * (u**10 + 11 * u**5 * v**5 - v**10)


def T(u, v):
    return (
        u**30
        + 522 * (u**25 * v**5 - u**5 * v**25)
        - 10005 * (u**20 * v**10 + u**10 * v**20)
        + v**30
    )


# Example function
def g(z):
    # return z**2  # try also 1/z, np.exp(z), etc.
    # return -(f(z, 1) ** 5) / T(z, 1) ** 2
    return -(H(z, 1) ** 3) / 12**3 / f(z, 1) ** 5


# Brightness (value) mapping for |f(z)|
# You can tweak "mag_scale" to compress/expand high magnitudes
def value_from_magnitude(r: np.ndarray, mag_scale: float = 1.0) -> np.ndarray:
    # Smoothly increases with |f|, approaches 1 as |f|→∞
    # mag_scale > 1.0 makes the map “brighter” for a given magnitude
    r_scaled = mag_scale * r
    return 1.0 - 1.0 / (1.0 + r_scaled)


# -----------------------------
# 2) Sphere parameterization
# -----------------------------
n_theta, n_phi = 200, 400
theta = np.linspace(0, np.pi, n_theta)
phi = np.linspace(0, 2 * np.pi, n_phi)
Theta, Phi = np.meshgrid(theta, phi)

X = np.sin(Theta) * np.cos(Phi)
Y = np.sin(Theta) * np.sin(Phi)
Z = np.cos(Theta)

# -----------------------------
# 3) Inverse stereographic projection (sphere → complex plane)
#    North pole maps to infinity
# -----------------------------
den = 1.0 - Z
Zc = np.full_like(X, np.inf, dtype=complex)  # default: ∞ at/near north pole
mask = den > 1e-10
Zc[mask] = (X[mask] + 1j * Y[mask]) / den[mask]

# -----------------------------
# 4) Evaluate function & make domain-coloring colors (HSV)
# -----------------------------
F = g(Zc)

# Hue from argument in [0,1]
hue = (np.angle(F) + np.pi) / (2.0 * np.pi)

# Value (brightness) from magnitude
mag = np.abs(F)
val = value_from_magnitude(mag, mag_scale=1.0)

# Saturation fixed at 1
sat = np.ones_like(val)

# Handle NaNs/Infs gracefully: set to max brightness, arbitrary hue
bad = ~np.isfinite(mag)
hue[bad] = 0.0  # red, but pick what you prefer
val[bad] = 1.0
sat[bad] = 1.0

RGB = hsv_to_rgb(np.stack((hue, sat, val), axis=-1))

# -----------------------------
# 5) Plot the sphere with domain coloring
# -----------------------------
fig = plt.figure(figsize=(9.5, 9))
ax = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(X, Y, Z, facecolors=RGB, linewidth=0, antialiased=False)

ax.set_box_aspect([1, 1, 1])
ax.set_axis_off()
ax.set_title("Domain coloring on the Riemann sphere", pad=12)

# -----------------------------
# 6) Add an Argument (Phase) Color Wheel Legend
# -----------------------------
# Build a unit-disk image colored by angle (hue). Brightness = 1 inside the disk, 0 outside.
res = 360
xx = np.linspace(-1, 1, res)
yy = np.linspace(-1, 1, res)
XX, YY = np.meshgrid(xx, yy)
ZZ = XX + 1j * YY
R = np.hypot(XX, YY)

wheel_hue = (np.angle(ZZ) + np.pi) / (2 * np.pi)
wheel_sat = np.ones_like(R)
wheel_val = (R <= 1).astype(float)

wheel_rgb = hsv_to_rgb(np.stack((wheel_hue, wheel_sat, wheel_val), axis=-1))

# Inset axes for the wheel (left lower corner); tweak [x0, y0, width, height] in figure coords
ax_wheel = fig.add_axes([0.07, 0.08, 0.22, 0.22])
ax_wheel.imshow(wheel_rgb, extent=(-1, 1, -1, 1), origin="lower")
# draw a circle boundary
theta_c = np.linspace(0, 2 * np.pi, 361)
ax_wheel.plot(np.cos(theta_c), np.sin(theta_c), linewidth=1)
ax_wheel.set_aspect("equal", "box")
ax_wheel.set_xticks([])
ax_wheel.set_yticks([])
ax_wheel.set_title("Argument (phase)")

# -----------------------------
# 7) Add a Magnitude → Brightness Bar
# -----------------------------
# Choose a representative magnitude range you care about
m_max = 5.0
m_samples = np.linspace(0, m_max, 400)
# Use a constant hue (e.g., blue) to visualize brightness only
bar_hue = np.full_like(m_samples, 0.66)  # 0.66 ~ blue
bar_sat = np.ones_like(m_samples)
bar_val = value_from_magnitude(m_samples, mag_scale=1.0)
bar_rgb = hsv_to_rgb(np.stack((bar_hue, bar_sat, bar_val), axis=-1))

ax_bar = fig.add_axes(
    [0.35, 0.13, 0.40, 0.05]
)  # [x0, y0, width, height] in figure coords
ax_bar.imshow(
    bar_rgb[np.newaxis, :, :], aspect="auto", extent=(0, m_max, 0, 1), origin="lower"
)
ax_bar.set_yticks([])
ax_bar.set_xlabel(r"$|f(z)|$")
ax_bar.set_title("Brightness vs. magnitude")

plt.show()
