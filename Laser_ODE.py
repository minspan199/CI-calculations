# -*- coding: utf-8 -*-
"""
Created on Wed May  9 14:30:08 2018

@author: Mingsen
"""
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp as ODE

Nd = 16
t1 = 1
t2 = 3
g = 0.5
h = -0.5493
N0 = 1e24
h = 0
tau_p = 0.02
pA = g * tau_p
pB = g * tau_p
dB = 0
tau_s = 2000 * tau_p
alpha = 3


def func(t, y):
    y1 = y0[0:Nd + 1]
    y2 = y0[Nd + 1:2 * Nd + 1]
    y3 = y0[2 * Nd + 1:3 * Nd + 2]
    dy1 = ((1 - 1j * alpha) * y3 * y1).reshape((1, -1)) - 1j * t1 * tau_p * np.exp(h) * np.append(y2, [
        0]) - 1j * t2 * tau_p * np.exp(-h) * np.append([0], y2)
    dy2 = -(pB - 1j * dB * tau_p) * y2 - 1j * t1 * tau_p * np.exp(-h) * y1[0:Nd] - 1j * t2 * tau_p * np.exp(h) * y1[
                                                                                                                 1:Nd + 1]
    dy3 = pA * tau_p / tau_s - tau_p / tau_s * y3 - tau_p / tau_s * (1 + 2 * y3) * np.abs(y1) * np.abs(y1)
    a = np.append(dy1.T, dy2)

    return np.append(a, dy3)


y0 = np.random.rand(3 * Nd + 2)
sol = ODE(func, [0, 5e4], y0)
