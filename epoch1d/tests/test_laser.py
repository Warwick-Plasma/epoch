#!/usr/bin/env python3

# Stephan Kuschel, 150818

import numpy as np
import sdf
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
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
    ax.set_title('{:2.1f} fs'.format(data['Header']['time']*1e15))

def plotevolution(key):
    sdffiles = ['{:04d}.sdf'.format(i) for i in range(8)]
    fig, axarr = plt.subplots(2,4, figsize=(16,9))
    for (sdffile, ax) in zip(sdffiles, np.ravel(axarr)):
        plotdump(sdffile, key, ax)
    fig.suptitle(key)
    fig.savefig(key.replace('/','_') + '.png', dpi=160)

def createplots():
    #showdatafields('0000.sdf')
    plotevolution('Electric Field/Ey')


class test_laser(SimTest):

    def totaleyassert(self, dump, val):
        data = sdf.read(dump, dict=True)
        ey = data['Electric Field/Ey'].data
        self.assertAlmostEqual(np.sum(ey), val)

    def test_createplots(self):
        createplots()

    def test_Eey0000(self):
        self.totaleyassert('0000.sdf', 0.0)

    def test_Eey0001(self):
        self.totaleyassert('0001.sdf', 206954574900.52179)

    def test_Eey0003 (self):
        self.totaleyassert('0003.sdf', 121431247301.71188)

    def test_Eey0007(self):
        self.totaleyassert('0007.sdf', -9490905.2679244578)

if __name__=='__main__':
    unittest.main()


