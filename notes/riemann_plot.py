import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, colors
from matplotlib.colors import hsv_to_rgb


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


def value_from_magnitude(magnitudes, scale=100.0):
    return 1 - np.exp(-scale * magnitudes)


# Cube-sphere parameterization
def cube_face_to_sphere(x_face, y_face, z_face):
    norm = np.sqrt(x_face**2 + y_face**2 + z_face**2)
    return x_face / norm, y_face / norm, z_face / norm


def make_cube_sphere_faces(n_face):
    grid = np.linspace(-1.0, 1.0, n_face)
    U, V = np.meshgrid(grid, grid)
    return [
        cube_face_to_sphere(np.ones_like(U), U, V),
        cube_face_to_sphere(-np.ones_like(U), U, V),
        cube_face_to_sphere(U, np.ones_like(U), V),
        cube_face_to_sphere(U, -np.ones_like(U), V),
        cube_face_to_sphere(U, V, np.ones_like(U)),
        cube_face_to_sphere(U, V, -np.ones_like(U)),
    ]


def stereographic_from_sphere(x_coord, y_coord, z_coord):
    denom = 1.0 - z_coord
    z_complex = np.full_like(x_coord, np.inf, dtype=complex)
    mask = denom > 1e-8
    z_complex[mask] = (x_coord[mask] + 1j * y_coord[mask]) / denom[mask]
    return z_complex


def colors_for_face(x_coord, y_coord, z_coord):
    z_complex = stereographic_from_sphere(x_coord, y_coord, z_coord)
    f_values = g(z_complex)
    magnitudes = np.abs(f_values)
    hue = (np.angle(f_values) + np.pi) / (2 * np.pi)
    val = value_from_magnitude(magnitudes)
    sat = np.ones_like(val)
    return hsv_to_rgb(np.stack((hue, sat, val), axis=-1))


n_face = 100
faces = make_cube_sphere_faces(n_face)

# Plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_axes([0.05, 0.08, 0.68, 0.84], projection="3d")
for X, Y, Z in faces:
    RGB = colors_for_face(X, Y, Z)
    ax.plot_surface(
        X,
        Y,
        Z,
        facecolors=RGB,
        linewidth=0,
        rstride=1,
        cstride=1,
        antialiased=False,
        shade=False,
    )
ax.set_box_aspect([1, 1, 1])
ax.set_title(r"Domain colouring on the Riemann sphere using $I(z)$")
ax.set_axis_off()

phase_cax = fig.add_axes([0.80, 0.18, 0.025, 0.64])
value_cax = fig.add_axes([0.90, 0.18, 0.025, 0.64])

phase_norm = colors.Normalize(vmin=-np.pi, vmax=np.pi)
phase_mappable = cm.ScalarMappable(norm=phase_norm, cmap=cm.hsv)
phase_mappable.set_array([])
cbar = fig.colorbar(phase_mappable, cax=phase_cax)
cbar.set_label(r"Hue ($\arg I(z)$)")
cbar.set_ticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
cbar.set_ticklabels([r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"])
cbar.ax.yaxis.labelpad = 10

value_max = 0.05
value_norm = colors.Normalize(vmin=0.0, vmax=value_max)
gray_values = np.linspace(0.0, 1.0, 256)
value_cmap = colors.ListedColormap(
    np.stack((gray_values, gray_values, gray_values), axis=-1)
)
value_mappable = cm.ScalarMappable(norm=value_norm, cmap=value_cmap)
value_mappable.set_array([])
value_cbar = fig.colorbar(value_mappable, cax=value_cax)
value_cbar.set_label(r"Brightness ($1-e^{-100|I(z)|}$)")
value_cbar.set_ticks([0.0, 0.01, 0.02, 0.03, 0.04, 0.05])
value_cbar.ax.yaxis.labelpad = 12

plt.show()
