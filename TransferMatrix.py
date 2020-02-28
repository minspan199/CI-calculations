# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 13:40:40 2018

@author: Mingsen
"""

# syms Km
import numpy as np

N = 20
ne = 3.3
L = 1
k = 4 * np.pi / L

Space = L / N
x = np.linspace(-0.5 * L, 0.5 * L, N + 1)
Ksi = (2 * ne * x / L) ** 2 - 1j * (2 * ne / k / L)
n = np.sqrt(Ksi)
Ms = np.identity(2)
D0 = np.mat([[1, 1], [ne * k, -ne * k]]).astype(complex)
Dn = np.mat([[1, 1], [ne * k, -ne * k]]).astype(complex)  # converting to complex 128

real = 200
imag = 600
Real = np.linspace(0, real - 1, real).astype(int)
Imaginary = np.linspace(0, imag - 1, imag).astype(int)
KT = np.linspace(0, N, N + 1).astype(int)
Z = np.zeros([real, imag])
for X in Real:
    for Y in Imaginary:
        Ms = np.identity(2)
        Km = (0.01 * (X + 1) - 0.001j * (Y - 300)) * k
        for kt in KT:
            Dv = np.mat([[1, 1], [n[kt] * Km, -n[kt] * Km]])
            Pv = np.mat([[np.exp(1j * n[kt] * Km * Space), 0], [0, np.exp(-1j * n[kt] * Km * Space)]])
            Mv = Dv * Pv * (Dv ** (-1))
            Ms = Ms * Mv
        Ms = D0 ** (-1) * Ms * Dn
        Z[X][Y] = np.abs(Ms[1, 1])
    print(X)
