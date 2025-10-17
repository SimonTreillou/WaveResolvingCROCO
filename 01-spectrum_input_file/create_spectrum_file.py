import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d

Nfrq=50*31
mode="analytical"
mode="DUCK94"
visual_check=True

#%%
if mode=="analytical":
    Tp=7.0
    fp=1/Tp
    sigma=0.01
    fmin=0.2*fp #maybe get inspired by SWASH?
    fmax=3.0*fp
    freq=np.linspace(fmin,fmax,Nfrq)
    S=np.exp(-(freq-fp)**2/sigma**2) # gaussian
    S=np.exp(-(freq-fp)**2/sigma**2) + 0.5 * np.exp(-(freq-1/15)**2/sigma**2) # double gaussian
elif mode=="DUCK94":
    day = '12'
    hh = '10'
    sens = 87
    df = pd.read_csv('/Users/simon/Data/Duck94_Elgar/Duck94_spectra/10' + day + hh + '00.ppp' + str(sens), header=None,
                     skiprows=5, sep=r'\s+')
    f1 = 0.01;f2 = 0.3
    if1 = np.argmin(np.abs(df[1] - f1))
    if2 = np.argmin(np.abs(df[1] - f2))
    f = interp1d(df[1], df[3]/1e4, kind='linear')
    freq = np.linspace(f1,f2,Nfrq)
    S = f(freq)*0.2
    df=freq[1]-freq[0]
    Hs=4*np.sqrt(np.sum(S*df))

if visual_check:
    plt.loglog(freq,S)
    plt.xlabel('frequency [Hz]')
    plt.ylabel('S')
    plt.show()

#%% Write
with open("./01-spectrum_input_file/FLASH_RIP/input_spectrum.txt", "w") as txt_file:
    for i in range(Nfrq):
        txt_file.write(str(freq[i])+" "+str(S[i])+"\n") # works with any number of elements in a line