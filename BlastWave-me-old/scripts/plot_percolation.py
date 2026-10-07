import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("../output/events/percolation_density.dat")
x, y, n = d.T
nx, ny = 100, 100
X = x.reshape(ny, nx)
Y = y.reshape(ny, nx)
N = n.reshape(ny, nx)

plt.figure(figsize=(6,5))
plt.pcolormesh(X, Y, N, shading='auto', cmap='hot')
plt.colorbar(label='string coverage')
plt.xlabel('x [fm]')
plt.ylabel('y [fm]')
plt.title('Percolation density map')
plt.tight_layout()
plt.savefig("../output/figures/phase5_percolation.png", dpi=140)
print("saved ../output/figures/phase5_percolation.png")
