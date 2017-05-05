#!/usr/bin/env python

# Copyright (C) 2016 Stephan Kuschel <Stephan.Kuschel@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import sdf
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import unittest
import platform
from . import SimTest


def showdatafields(sdffile):
    data = sdf.read(sdffile, dict=True)
    print('Data fields present in "' + sdffile + '":')
    print(data.keys())

def plotdump(sdffile, key, ax):
    data = sdf.read(sdffile, dict=True)
    axis = data[key].grid_mid.data[0]*1e6
    array = data[key].data
    ax.plot(axis, array)
    ax.set_title('{:2.1f} fs'.format(data['Header']['time']*1e15))
    ax.set_xlabel('x [µm]')

def plotevolution(key):
    sdffiles = ['{:04d}.sdf'.format(i) for i in range(8)]
    fig, axarr = plt.subplots(2,4, figsize=(16,9))
    for (sdffile, ax) in zip(sdffiles, np.ravel(axarr)):
        plotdump(sdffile, key, ax)
    fig.suptitle(key)
    fig.savefig(key.replace('/','_') + '.png', dpi=160)

def createplots():
    if platform.system() == 'Darwin':
        print('macosx backend')
        plt.switch_backend('macosx')
    #showdatafields('0000.sdf')
    plotevolution('Electric Field/Ey')


class test_maxwell_solvers(SimTest):

    def test_createplots(self):
        os.chdir('yee')
        createplots()

        os.chdir('../lehe_x')
        createplots()


if __name__=='__main__':
    unittest.main()


