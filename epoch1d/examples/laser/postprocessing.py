#!/usr/bin/env python

# Stephan Kuschel, 150818

import numpy as np
import sdf
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

def showdatafields(sdffile):
    data = sdf.read(sdffile, dict=True)
    print('Data fields present in "' + sdffile + '":')
    print(data.keys())

def plotdump(sdffile, key, ax):
    data = sdf.read(sdffile, dict=True)
    array = data[key].data
    ax.plot(array)
    ax.set_title('{:2.1f} fs'.format(data['Header']['time']*1e15))

def plotevolution(key):
    sdffiles = ['{:04d}.sdf'.format(i) for i in range(8)]
    fig, axarr = plt.subplots(2,4, figsize=(16,9))
    for (sdffile, ax) in zip(sdffiles, np.ravel(axarr)):
        plotdump(sdffile, key, ax)
    fig.suptitle(key)
    fig.savefig(key.replace('/','_') + '.png', dpi=160)



if __name__=='__main__':
    showdatafields('0000.sdf')
    plotevolution('Electric Field/Ey')
