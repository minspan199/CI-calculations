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

# clc
# clear all
# close all
# global Nd t1 t2 g h N0 tau_p pA pB dB tau_s alpha
# Nd = 16 t1 = 1 t2 = 3 g = 0.5 h = -0.5493 N0 = 1e24h = 0
# tau_p = 0.02 pA = g*tau_p pB = g*tau_p dB = 0 tau_s = 2000*tau_p alpha = 3
# T_span = [0 50000] y0 = rand(3*Nd + 2,1)
# [T,Y] = ode45(@rate_equation,T_span,y0)
#
#
# figure
# N_span = length(T)
# for x0 = 1:1:Nd + 1
#    n0 = 1:10:N_span
#    Yi = angle(Y(n0,x0))
#    p = plot3(T(n0),x0*ones(length(n0),1),Yi)
#    hold on
# end
# set(gcf, 'Position', [00, 00, 350, 300])
# zlim([-3.5,3.5])
# xlabel('Time')
# ylabel('Site number')
# zlabel('Phase')
# % subplot(2,1,2)
#
#
# figure
# % subplot(2,1,1)
# plot(T,abs(Y(:,1:Nd + 1)))
# ylim([-0.1 0.2])
# title('Density of carriers in the upper laser level') % carriers density in high laser level
# xlabel('Time')
# ylabel('Density')
# set(gcf, 'Position', [00, 00, 350, 300])
# set(gca,'FontSize', 14) % Font Size
#
# figure
# Hf = surf(abs(Y(:,1:Nd + 1)))
# set(gcf, 'Position', [00, 00, 350, 300])
# set(Hf,'LineStyle','none')
# zlim([0,0.5])
# ylabel('Time')
# xlabel('Site number')
# zlabel('Carrier densities')
#
# colormap(jet)
#
#
# figure
# N_span = length(T)
# for x0 = 1:1:Nd + 1
#    n0 = 1:10:N_span
#    Yi = abs(Y(n0,x0))
#    plot3(T(n0),x0*ones(length(n0),1),Yi,'LineWidth',2)
#    hold on
# end
# set(gcf, 'Position', [00, 00, 350, 300])
# zlim([0,0.5])
# xlabel('Time')
# ylabel('Site number')
# zlabel('Carrier densities')
# % subplot(2,1,2)
#
#
# figure
# plot(T,angle(Y(:,1:Nd + 1)))
# title('Density of the photons in the active medium') % photons density in activer region And this is the function.
# set(gca,'FontSize', 14) % Font Size
# set(gcf, 'Position', [00, 00, 350, 300])
# xlabel('Time')
# ylabel('Phase')
#
# figure
# subplot(2,1,1)
# plot(T,abs(Y(:,Nd + 2:2*Nd + 1)))
# ylim([-1 2])
# title('Densità dei portatori nel livello laser superiore') % carriers density in high laser level
# subplot(2,1,2)
# plot(T,angle(Y(:,Nd + 2:2*Nd + 1)))
# title('Densità dei fotoni nel mezzo attivo') % photons density in activer region And this is the function.
#
#
# function dy = rate_equation(t,y)
#
# global Nd t1 t2 g h N0 tau_p pA pB dB tau_s alpha
# y1 = y(1:Nd + 1)
# y2 = y(Nd + 2:2*Nd + 1)
# y3 = y(2*Nd + 2:3*Nd + 2)
# dy1 = (1 - 1i*alpha)*y3.*y1 - 1i*t1*tau_p*exp(h)*[y20] - 1i*t2*tau_p*exp(-h)*[0y2]
# dy2 = -(pB - 1i*dB*tau_p)*y2 - 1i*t1*tau_p*exp(-h)*y1(1:Nd) - 1i*t2*tau_p*exp(h)*y1(2:Nd + 1)
# dy3 = pA*tau_p/tau_s - tau_p/tau_s*y3 - tau_p/tau_s*(1 + 2*y3).*abs(y1).*abs(y1)
# dy = [dy1dy2dy3]
# end
#
# % this is a very simple (dimensionless) and efficient version of the rate equations (note you % can add another equation, dy(3)(not coupled), that gives you the phase(t) of your laser!
