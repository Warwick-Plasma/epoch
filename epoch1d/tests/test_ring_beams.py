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
import unittest
import platform
#from . import SimTest


def showdatafields(sdffile):
    data = sdf.read(sdffile, dict=True)
    print('Data fields present in "' + sdffile + '":')
    print(data.keys())


def save_plot(key):
    sdffile = '0000.sdf'
    fig = plt.figure()
    ax = plt.gca()
    data = sdf.read(sdffile, dict=True)
    x = data['Grid/' + key].data[0]
    y = data['dist_fn/' + key].data
    ax.plot(x, y)
    ax.set_title(key + ' ' + r'{:2.1f} ms'.format(data['Header']['time']*1e3))
    plt.tight_layout()
    plt.grid(b=True, which='both')
    fig.savefig(key.replace('/', '_') + '.png', dpi=160)


def save_contourf(key):
    sdffile = '0000.sdf'
    fig = plt.figure()
    ax = plt.gca()
    data = sdf.read(sdffile, dict=True)
    x = data['Grid/' + key].data[0]
    y = data['Grid/' + key].data[0]
    z = data['dist_fn/' + key].data[:, :]
    ax.contourf(x, y, z)
    ax.set_title(key + ' ' + r'{:2.1f} ms'.format(data['Header']['time']*1e3))
    plt.grid()
    fig.savefig(key.replace('/', '_') + '.png', dpi=160)


def createplots():
    if platform.system() == 'Darwin':
        print('macosx backend')
        plt.switch_backend('macosx')

    #showdatafields('0000.sdf')
    save_plot('px/ion_ring_beam')
    save_plot('py/ion_ring_beam')
    save_plot('pz/ion_ring_beam')
    save_contourf('px_py/ion_ring_beam')
    save_contourf('px_pz/ion_ring_beam')
    save_contourf('py_pz/ion_ring_beam')


class test_ring_beam(unittest.TestCase):#SimTest):

    def test_createplots(self):
        createplots()


if __name__ == '__main__':
    unittest.main()
