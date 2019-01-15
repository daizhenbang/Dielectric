import numpy as np
import matplotlib.pyplot as plt
import math as ma
from scipy.fftpack import fft

file = open('all_dielectric.txt','r');
#file = open('dielectric_cubic.txt','r');
lines = file.readlines();
file.close();
temperature = 77;

numSnapshot = int(len(lines)/3);
dielectric = np.zeros((numSnapshot,9));
count = 0;
for i in range(numSnapshot):
    for j in range(3):
        temp = lines[count].split();
        dielectric[i,j*3:(j+1)*3] = temp;
        count += 1;


lattice = np.zeros((3,3));
'''Orthorhombic'''
lattice[0][0] = 8.58;
lattice[1][1] = 8.87;
lattice[2][2] = 12.63;
'''Cubic'''
#lattice[0][0] = 8.925809001399999;
#lattice[1][1] = 8.925809001399999;
#lattice[2][2] = 12.623000144000001;
volume = abs(np.linalg.det(lattice));

polarizability = np.zeros((numSnapshot,9));
for i in range(numSnapshot):
    for j in range(9):
        if j == 0 or j == 4 or j == 8:
            polarizability[i,j] = (1-dielectric[i,j])*volume/4/ma.pi;
        else:
            polarizability[i,j] = (0-dielectric[i,j])*volume/4/ma.pi;

avePolariza = np.mean(polarizability,axis=0);
for i in range(numSnapshot):
    polarizability[i,:] = polarizability[i,:] - avePolariza;

'''Calculate the autocorrelation function'''
step = 1;
periodHighLimit = 200;
numAuto = int(periodHighLimit/step);
autoPolar0 = np.zeros((numAuto,9));


for i in range(numAuto):
    window = (i+1)*step;
    for j in range(9):
        for k in range(numSnapshot-window):
            autoPolar0[i,j] += polarizability[k,j]*polarizability[k+window,j];
        autoPolar0[i,j] = autoPolar0[i,j]/(numSnapshot-step);

'''Zero padding'''
numPadding = 200;
zeros0 = np.zeros((numPadding,9));
autoPolar = np.concatenate((autoPolar0,zeros0),axis=0);
numAuto += numPadding;
'''End of calculating the autocorrelation function'''

#plt.plot(autoPolar[:,0])
'''Output the autocorrelation. C, N, H, Pb, I'''
powerH = fft(autoPolar[:,0]+autoPolar[:,4]+autoPolar[:,8]);
fmax = 1/(100*10**(-15))*3.335641*10**(-11);
fmin = fmax/numAuto;

freq = np.arange(fmin,fmax+fmin,fmin);
kb = 8.6173303*10**(-5);
for i in range(numAuto):
    powerH[i] = powerH[i]*freq[i]*(1+1/(np.exp(freq[i]*0.0001239/temperature/kb)-1));

raman = np.abs(powerH);
for i in range(numAuto):
    print(str(freq[i])+"  "+str(raman[i]))

'''Do the broadening'''
#xgrid = np.arange(1,160,0.01);
#broadened = np.zeros((len(xgrid),));
#gauss = np.zeros((numAuto,len(xgrid)));
#sigma = 0.7;
#for i in range(numAuto):
#    for j in range(len(xgrid)):
#        gauss[i,j] = raman[i]*1/(2*ma.pi)**(0.5)/sigma*ma.exp(-(xgrid[j] - freq[i])**2/2/sigma**2);
#
#for i in range(numAuto):
#    broadened = broadened + gauss[i,:];
#
#plt.plot(xgrid,broadened)
'''End of broadening'''
#plt.plot(freq,raman)
#plt.xlim([5,fmax/2]);
#plt.ylim([0,8000000]);

