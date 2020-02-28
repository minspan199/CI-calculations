# -*- coding: utf-8 -*-
"""
Created on Sat May  5 13:17:47 2018

@author: Mingsen
"""

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
"import tensorflow as tf"

#def clear_all():
#    """Clears all the variables from the workspace of the spyder application."""
#    gl = globals().copy()
#    for var in gl:
#        if var[0] == '_': continue
#        if 'func' in str(globals()[var]): continue
#        if 'module' in str(globals()[var]): continue
#        del globals()[var]
#    if __name__ == "__untitled0__":
#        clear_all()
#clear_all()

#Parameter set
Nd = 16; t1 = 1; t2 = 2; g = 0.5; h = -0.5493; h = -0.3466;

#Hamiltonian assembling
H0 = np.diag(np.append(np.tile([1j*g,-1j*g],Nd),1j*g));
H1 = np.diag(np.tile([t1*np.exp(h),t2*np.exp(h)],Nd),1);
H2 = np.diag(np.tile([t1*np.exp(-h),t2*np.exp(-h)],Nd),-1);
H = H0 + H1 + H2; 
[w, v] = LA.eig(H)
L = np.linspace(0, 2*Nd, num = 2*Nd + 1).astype(int);
L0 = np.linspace(0, 2*Nd + 2, num = 2*Nd + 3).astype(int);
for i in L:
    if np.abs(w[i].real) < 0.0001:
        S_Mode_Index = [i];
state_amplitude = v[:,S_Mode_Index];
plt.figure(num = None, figsize = (4, 3), dpi = 120, facecolor = 'w', edgecolor = 'k')
index = np.arange(len(state_amplitude));
y_pos = []
for i in L:
    y_pos.append(np.real(state_amplitude)[i][0]);

plt.bar(index + 1, y_pos)
plt.plot(L0, np.zeros(2*Nd + 3))
plt.xlim(0, 2*Nd + 2)
plt.title('Eigenstates of the zero mode')
plt.xlabel('Site number')
plt.ylabel('Amplitudes')
plt.xlabel('Site number')
plt.xticks([1, Nd/2 + 1, Nd + 1, Nd*1.5 + 1, 2*Nd +1])

#H = diag([repmat([1i*g -1i*g],[1 Nd]) 1i*g]) +...
#    diag(repmat([t1*exp(h) t2*exp(h)],[1 Nd]),1) +...
#    diag(repmat([t1*exp(-h) t2*exp(-h)],[1 Nd]),-1);  % Assembling of Hamiltonian matrix
#[V, A] = eig(H);
#lam = diag(A);
#[~, idx] = sort(real(lam));
#lam1 = lam(idx);
#figure
#plot(real(lam1),'b*')
#hold on
#plot(imag(lam1),'r*')
#set(gcf, 'Position', [00, 00, 350, 300])
#axis([0 2*Nd + 2 -4 4])
#set(gca,'FontSize', 14) % Font Size
#figure
#bar(V(:,17))
#set(gcf, 'Position', [00, 00, 350, 300])
#set(gca,'FontSize', 14) % Font Size
#axis([0 2*Nd + 2 -1 1])

#H = np.array([[1+1j,2+1j, 3+1j],[3+2j, 2+2j, 1+1j],[3+8j, 4+2j, 6-1j]])
#[a, b] = tf.eig(H)
