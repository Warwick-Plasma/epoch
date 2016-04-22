#!/usr/bin/env python

# Stephan Kuschel, 150818

import numpy as np
import sdf
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import unittest
from . import SimTest


def showdatafields(sdffile):
    data = sdf.read(sdffile, dict=True)
    print('Data fields present in "' + sdffile + '":')
    print(data.keys())

def plotdump(sdffile, key, ax):
    data = sdf.read(sdffile, dict=True)
    array = data[key].data
    ax.plot(array)
    ax.set_title(r'{:2.1f} ms'.format(data['Header']['time']*1e3))

def plotdump2d(sdffile, key, ax):
    data = sdf.read(sdffile, dict=True)
    array = data[key].data[:,:,0]
    ax.imshow(array)
    ax.set_title(r'{:2.1f} ms'.format(data['Header']['time']*1e3))

def plotevolution(key):
    print('plotting: ' + key)
    sdffiles = ['{:04d}.sdf'.format(i) for i in range(0,200,14)]
    fig, axarr = plt.subplots(3,5, figsize=(20,11))
    for (sdffile, ax) in zip(sdffiles, np.ravel(axarr)):
        plotdump(sdffile, key, ax)
    fig.suptitle(key)
    fig.savefig(key.replace('/','_') + '.png', dpi=160)

def plot2devolution(key):
    print('plotting: ' + key)
    sdffiles = ['{:04d}.sdf'.format(i) for i in range(0,200,14)]
    fig, axarr = plt.subplots(3,5, figsize=(20,11))
    for (sdffile, ax) in zip(sdffiles, np.ravel(axarr)):
        plotdump2d(sdffile, key, ax)
    fig.suptitle(key)
    fig.savefig(key.replace('/','_') + '.png', dpi=160)

def createplots():
    #showdatafields('0000.sdf')
    plotevolution('Electric Field/Ex')
    plotevolution('Derived/Number_Density/electrons')
    plotevolution('Derived/Number_Density/protons')
    plotevolution('Derived/Charge_Density')
    plot2devolution('dist_fn/x_px/electrons')
    plot2devolution('dist_fn/x_px/protons')


class test_landau(SimTest):

    def test_creatplots(self):
        createplots()


if __name__=='__main__':
    unittest.main()


