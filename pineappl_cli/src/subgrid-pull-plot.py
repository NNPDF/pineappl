#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from math import fabs, log10
from scipy.interpolate import griddata

x1 = np.array([{x1}])  # noqa: F821
x2 = np.array([{x2}])  # noqa: F821
z = np.array([{z}])  # noqa: F821
x = 0.5 * np.log(x1 / x2)
y = np.sqrt(x1 * x2)

nrap = 50
nmass = 50

sym_min = -max(fabs(np.min(x)), fabs(np.max(x)))
sym_max = max(fabs(np.min(x)), fabs(np.max(x)))

xi = np.linspace(sym_min, sym_max, (nrap // 2) * 2 + 1)
yi = np.logspace(log10(np.min(y)), log10(np.max(y)), nmass)
zi = griddata(
    (x, y), z, (xi[None, :], yi[:, None]), method="linear", rescale=True
)

# print(xi.shape)
# print(yi.shape)
# print(zi.shape)

# mask impossible kinematic values
for iy, ix in np.ndindex(zi.shape):
    # print(ix, iy)
    x1v = yi[iy] * np.exp(xi[ix])
    x2v = yi[iy] / np.exp(xi[ix])

    # print('y = {{}} m/s = {{}} -> x1 = {{}} x2 = {{}}'.format(xi[ix], yi[iy], x1v, x2v))

    if x1v > 1.0 or x2v > 1.0:
        zi[iy, ix] = np.nan

figure, axes = plt.subplots(1, 2, constrained_layout=True)
figure.set_size_inches(10, 5)

mesh = axes[0].pcolormesh(xi, yi, zi, shading="nearest", linewidth=0, snap=True)
axes[0].scatter(x, y, marker="*", s=5)
axes[0].set_yscale("log")
axes[0].set_xlabel(r"$y = 1/2 \log (x_1/x_2)$")
axes[0].set_ylabel(r"$M/\sqrt{{s}} = \sqrt{{x_1 x_2}}$")
# axes[0].set_aspect('equal', share=True)

x1i = np.logspace(log10(np.min(x1)), log10(np.max(x1)), 50)
x2i = np.logspace(log10(np.min(x2)), log10(np.max(x2)), 50)
z12i = griddata(
    (x1, x2), z, (x1i[None, :], x2i[:, None]), method="linear", rescale=True
)

mesh = axes[1].pcolormesh(
    x1i, x2i, z12i, shading="nearest", linewidth=0, snap=True
)
axes[1].set_xscale("log")
axes[1].set_yscale("log")
axes[1].scatter(x1, x2, marker="*", s=5)
axes[1].set_aspect("equal", share=True)
axes[1].set_xlabel(r"$x_1$")
axes[1].set_ylabel(r"$x_2$")

figure.colorbar(mesh, ax=axes, extend="min")
figure.savefig("plot.pdf")
