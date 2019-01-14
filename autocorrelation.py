import numpy as np
import matplotlib.pyplot as plt
import math as ma
from scipy.fftpack import fft

file = open('all_dielectric.txt','r');
lines = file.readlines();
file.close();

numSnapshot = int(len(lines)/3);
dielectric = np.zeros((numSnapshot,9));
count = 0;
for i in range(numSnapshot):
    for j in range(3):
        temp = lines[count].split();
        dielectric[i,j*3:(j+1)*3] = temp;
        count += 1;


lattice = np.zeros((3,3));
lattice[0][0] = 8.58;
lattice[1][1] = 8.87;
lattice[2][2] = 12.63;
volume = abs(np.linalg.det(lattice));

polarizability = np.zeros((numSnapshot,9));
for i in range(numSnapshot):
    for j in range(9):
        if j == 0 or j == 4 or j == 8:
            polarizability[i,j] = (1-dielectric[i,j])*volume/4/ma.pi;
        else:
            polarizability[i,j] = (0-dielectric[i,j])*volume/4/ma.pi;


'''Calculate the autocorrelation function'''
step = 1;
periodHighLimit = 299;
numAuto = int(periodHighLimit/step);
autoPolar = np.zeros((numAuto,9));


for i in range(numAuto):
    window = (i+1)*step;
    for j in range(9):
        for k in range(numSnapshot-window):
            autoPolar[i,j] += polarizability[k,j]*polarizability[k+window,j];
        autoPolar[i,j] = autoPolar[i,j]/(numSnapshot-window);

'''End of calculating the autocorrelation function'''

'''Output the autocorrelation. C, N, H, Pb, I'''
powerH = fft(autoPolar[:,0]+autoPolar[:,4]+autoPolar[:,8]);
fmax = 1/(100*10**(-15))*3.335641*10**(-11);
fmin = fmax/numAuto;

freq = np.arange(0,fmax,fmin);
kb = 8.6173303*10**(-5);
for i in range(numAuto):
    powerH[i] = powerH[i]*freq[i]*(1+1/(np.exp(freq[i]*0.0001239/77/kb))-1);

plt.plot(freq,np.abs(powerH))
plt.xlim([5,fmax/2]);
plt.ylim([0,500000])
