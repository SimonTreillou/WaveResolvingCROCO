import numpy as np
import matplotlib.pyplot as plt

Nfrq=50*31
Tp=7.0
fp=1/Tp
sigma=0.01
fmin=0.2*fp #maybe get inspired by SWASH?
fmax=3.0*fp
freq=np.linspace(fmin,fmax,Nfrq)
S=np.exp(-(freq-fp)**2/sigma**2) # gaussian

S=np.exp(-(freq-fp)**2/sigma**2) + 0.5 * np.exp(-(freq-1/15)**2/sigma**2) # double gaussian

plt.plot(freq,S)
plt.show()

with open("./input_spectrum.txt", "w") as txt_file:
    for i in range(Nfrq):
        txt_file.write(str(freq[i])+" "+str(S[i])+"\n") # works with any number of elements in a line
