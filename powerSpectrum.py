import numpy as np
import matplotlib.pyplot as plt
import math as ma
from scipy.fftpack import fft

#def GaussianBroaden(data,sigma):
#    size = len(data);
#    broadened = np.zeros((size,1));
#    base = np.arange(0,size,1);
#    for i in range(size):
#        temp0 = [x - i for x in base];
#        temp = np.square(temp0)/2/sigma/sigma;
##        print(temp)
##        for j in range(size):
##            temp[i] = ma.exp(temp[i]);
#        broadened = broadened + data[i]*np.exp(temp);
#        
#    return broadened;
    

'''import the autocorelation'''
name = 'C';

#['I','Pb','C','N','H']:
mass = [126.90,207.2,12.0100,14.010,1.0080];
filename = 'auto'+name+'.txt'
file1 = open(filename,'r');
lines = file1.readlines();
file1.close;

length = len(lines);

autoH = np.zeros((length,));


for i in range(length):
    autoH[i] = float(mass[2])*float(lines[i]);
    
'''Sum up all the autocorrelation'''
autoH = np.zeros((length,));
mass = [126.90,207.2,12.0100,14.010,1.0080];
count = 0;
for name in ['I','Pb','C','N','H']:
    filename = 'auto'+name+'.txt'
    file1 = open(filename,'r');
    lines = file1.readlines();
    file1.close;

    length = len(lines);    
    for i in range(length):
        autoH[i] = autoH[i] + float(mass[count])*float(lines[i]);
    
    count += 1;
'''End of it'''

time = np.arange(0,length,1);
time = time*5;

powerH = fft(autoH);
fmax = 1/(5*10**(-15))*3.335641*10**(-11);
fmin = fmax/length;

freq = np.arange(0,fmax,fmin);
freq = freq;

fig = plt.figure(1,figsize = (20,15))
sub1 = fig.add_subplot(2,2,1);
sub1.plot(freq,np.abs(powerH))
plt.xlim([0,200])
plt.ylim([-0.001,0.04])
plt.xlabel(r'$cm^{-1}$')
#plt.title(name)
plt.title('All Atoms')

sub2 = fig.add_subplot(2,2,2);
sub2.plot(freq,np.abs(powerH));
plt.xlim([0,fmax/2])
plt.xlabel(r'$cm^{-1}$')

sub3 = fig.add_subplot(2,2,3);
sub3.plot(time,autoH);
plt.xlabel(r'$time/fs$')

picName= name+'.png';
picName = 'all.png'
fig.savefig(picName)