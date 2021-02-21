#!/usr/bin/env python
# By KS 01/21 (a first python attempt)
# Plot the cross sections from the specified channel file

import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

font = {'family' : 'Times',
        'size'   : 18}
plt.rc('font', **font)

if len(sys.argv) != 2 :
    print("Usage: ./xscn_plotter.py <channel file name> ")
    sys.exit()
else:
    channelfile = sys.argv[1]

#channelfile = "argon"
chanfilename = "../channels/channels_"+channelfile+".dat"
chans = np.genfromtxt(chanfilename,dtype={'names': ('cname', 'num', 'cpstate','flavor','ntargets'),'formats': ('U16', 'i', 'U12','U12','i')})


ich=0
while ich<chans.size:
    channame = chans['cname'][ich]
    xscnfilename = "../xscns/xs_"+str(channame)+".dat"
    print(xscnfilename)
    xscn = np.loadtxt(xscnfilename, comments='#', skiprows=0)
    energyingev = 10**xscn[:,0]
    energyinmev = energyingev*1000.
# Select the right column in the cross section file
    colnum=1

    if chans['flavor'][ich] == "e" and chans['cpstate'][ich] == "+":
        colnum=1
    if chans['flavor'][ich] == "m" and chans['cpstate'][ich] == "+":
        colnum=2
    if chans['flavor'][ich] == "t" and chans['cpstate'][ich] == "+":
        colnum=3
    if chans['flavor'][ich] == "e" and chans['cpstate'][ich] == "-":
        colnum=4
    if chans['flavor'][ich] == "m" and chans['cpstate'][ich] == "-":
        colnum=5
    if chans['flavor'][ich] == "t" and chans['cpstate'][ich] == "-":
        colnum=6
    xscnval = xscn[:,colnum]*energyingev
    plt.plot(energyinmev,xscnval,label=chans['cname'][ich])
    
    ich=ich+1
    
plt.yscale('log')
plt.xlabel("Neutrino energy (MeV)")
plt.ylabel("Cross section ($10^{-38}$ cm$^2$)")
plt.legend(loc='center right',prop={'size': 12})
#bottom, top = plt.ylim()
#plt.ylim(bottom, top)
left, right = plt.xlim()
plt.xlim(left, right*1.5)
plt.show()
