import math

import h5py
import matplotlib.pyplot as plt
import numpy as np

with h5py.File('spectrum-5.mat', 'r') as f:
    f.keys()
    print(list(f.keys()))
    Z_22 = f.get('Z_22')
    data = np.array(Z_22)
    data_inv = np.array(Z_22)

for ix, iy in np.ndindex(data.shape):
    data_inv[ix, iy] = math.log(abs(1 / data[ix, iy]))

fig, axes = plt.subplots(nrows=1, ncols=1)
spect = axes.imshow(data_inv, origin='lower', interpolation='nearest', aspect='auto', extent=[0.01, 1.9, -0.35, 0.35],
                    vmin=-4, vmax=7, cmap='jet')
# Add minorticks on the colorbar to make it easy to read the
# values off the colorbar.
cbar = fig.colorbar(spect, ax=axes, extend='both')
cbar.minorticks_on()
plt.title('eigenmodes')
plt.ylabel('Imaginary part')
plt.xlabel('Wavenumber')
plt.show()
