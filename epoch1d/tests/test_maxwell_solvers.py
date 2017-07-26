#!/usr/bin/env python

# Copyright (C) 2017 Alexander Blinne <A.Blinne@gsi.de>,
# Stephan Kuschel <Stephan.Kuschel@gmail.com>
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

import types
import os
import os.path as osp
import unittest
import platform

import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sdf

from . import SimTest


micron = 1e-6
femto = 1e-15
c = 2.99792458e8 # m/s

# check that these correspond to input deck!

nx = 240

x_min = -12 * micron
x_max = -x_min

dt_multiplier = 0.95

lambda_l = 0.5 * micron
x0 = -12.0 * micron # m
t0 = 8 * femto # s

# derived quantities from the above
k_l = 2*np.pi/lambda_l
dx = (x_max-x_min)/nx

dt = dt_multiplier * dx / c

vg_lehe = c*(1.0 + 2.0*(1.0-c*dt/dx)*(k_l*dx/2.0)**2)
vg_yee = c*np.cos(k_l*dx/2.0)/np.sqrt(1-(c*dt/dx*np.sin(k_l*dx/2.0))**2)


def xt(sdffile, key = 'Electric Field/Ey'):
    t = sdffile['Header']['time']
    xaxis = sdffile[key].grid_mid.data[0]
    data = sdffile[key].data
    b = np.sum(data**2)
    if b>0:
        x = np.sum(xaxis*data**2)/b
    else:
        x = None

    return t, x


class test_maxwell_solvers(SimTest):
    @classmethod
    def setUpClass(cls):
        super(test_maxwell_solvers, cls).setUpClass()
        sdffiles = zip([osp.join('yee','{:04d}.sdf'.format(i)) for i in range(8)],   [osp.join('lehe_x','{:04d}.sdf'.format(i)) for i in range(8)])

        dumps = cls.dumps = []
        for yee_dump, lehe_dump in sdffiles:
            dump = types.SimpleNamespace()
            dump.yee = sdf.read(yee_dump, dict=True)
            dump.lehe = sdf.read(lehe_dump, dict=True)
            dumps.append(dump)


    def test_createplot(self):
        if platform.system() == 'Darwin':
            print('macosx backend')
            plt.switch_backend('macosx')

        fig, axarr = plt.subplots(2,4, figsize=(16,9))

        key = 'Electric Field/Ey'

        for (dump, ax) in zip(self.dumps, np.ravel(axarr)):
            array_yee = dump.yee[key].data
            array_lehe = dump.lehe[key].data
            axis = dump.yee[key].grid_mid.data[0]*1e6

            ax.plot(axis, array_lehe, label='Lehe', linewidth=1)
            ax.plot(axis, array_yee, label='Yee', linewidth=1)
            ax.set_title('{:2.1f} fs'.format(dump.yee['Header']['time']*1e15))
            ax.set_xlabel(r'x [${\mu}\mathrm{m}$]')
            ax.legend()
        fig.suptitle(key)
        fig.savefig(key.replace('/','_') + '.png', dpi=320)

    def test_group_velocity(self):
        tx_yee = []
        tx_lehe = []
        for dump in self.dumps:
            t, x = xt(dump.yee)
            if t < 6e-15:
                continue
            tx_yee.append((t,x))
            tx_lehe.append(xt(dump.lehe))

        tx_yee = np.array(tx_yee)
        vg_yee_sim = np.polyfit(tx_yee[:,0], tx_yee[:,1], 1)[0]
        print(tx_yee)

        # for reference, right here, right now, the following line prints
        # 292363351.796 291329547.371 0.00353602604066
        print(vg_yee, vg_yee_sim, abs(vg_yee-vg_yee_sim)/vg_yee)

        tx_lehe = np.array(tx_lehe)
        vg_lehe_sim = np.polyfit(tx_lehe[:,0], tx_lehe[:,1], 1)[0]
        print(tx_lehe)

        # for reference, right here, right now, the following line prints
        # 311627789.85156083 310055314.605 0.00504600455477
        print(vg_lehe, vg_lehe_sim, abs(vg_lehe-vg_lehe_sim)/vg_lehe)

        assert np.isclose(vg_yee, vg_yee_sim, rtol=0.01) #
        assert np.isclose(vg_lehe, vg_lehe_sim, rtol=0.01)

if __name__=='__main__':
    unittest.main()


