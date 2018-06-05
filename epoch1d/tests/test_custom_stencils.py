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
c = 2.99792458e8  # m/s

# check that these correspond to input deck!

nx = 240

x_min = -12 * micron
x_max = -x_min

dt_multiplier = 0.95

lambda_l = 0.5 * micron
x0 = -12.0 * micron  # m
t0 = 8 * femto  # s

# derived quantities from the above
k_l = 2*np.pi/lambda_l
dx = (x_max-x_min)/nx

dt = dt_multiplier * dx / c

vg_lehe = c*(1.0 + 2.0*(1.0-c*dt/dx)*(k_l*dx/2.0)**2)
vg_yee = c*np.cos(k_l*dx/2.0)/np.sqrt(1-(c*dt/dx*np.sin(k_l*dx/2.0))**2)

# calculate_omega.py --dim 2 --symmetric --params 0.95 0.125
# -0.013364149548965119 --physical-dx 1e-7 --laser-wavelength 5e-7
# yields vg=1.0062495084969005 c
vg_opt = c * 1.0062495084969005


def xt(sdffile, key='Electric Field/Ey'):
    t = sdffile['Header']['time']
    xaxis = sdffile[key].grid_mid.data[0]
    data = sdffile[key].data
    b = np.sum(data**2)
    if b > 0 and t > 0:
        x = np.sum(xaxis*data**2)/b
    else:
        x = None

    return t, x


class test_custom_stencils(SimTest):
    solvers = ['lehe_x', 'lehe_custom', 'optimized']

    @classmethod
    def setUpClass(cls):
        super(test_custom_stencils, cls).setUpClass()

        dumps = cls.dumps = {}
        for solver in cls.solvers:
            l = dumps.setdefault(solver, [])
            for dump in [osp.join(solver, '{:04d}.sdf'.format(i))
                         for i in range(8)]:
                l.append(sdf.read(dump, dict=True))

    def test_createplot(self):
        if platform.system() == 'Darwin':
            print('macosx backend')
            plt.switch_backend('macosx')

        key = 'Electric Field/Ey'
        fig, axarr = plt.subplots(2, 4, figsize=(16, 9))

        for (dump_id, ax) in enumerate(np.ravel(axarr)):
            data_lehe_custom = self.dumps['lehe_custom'][dump_id]
            data_lehe_x = self.dumps['lehe_x'][dump_id]
            data_optimized = self.dumps['optimized'][dump_id]

            axis = data_optimized[key].grid_mid.data[0]*1e6
            array_optimized = data_optimized[key].data
            array_lehe_x = data_lehe_x[key].data
            array_lehe_custom = data_lehe_custom[key].data
            array_optimized = data_optimized[key].data
            ax.plot(axis, array_lehe_x, label='Lehe (Builtin)', linewidth=1)
            ax.plot(axis, array_lehe_custom, label='Lehe (Custom)',
                    linewidth=1, linestyle='dotted')
            ax.plot(axis, array_optimized, label='Optimized', linewidth=1)
            ax.set_title('{:2.1f} fs'.format(
                         data_optimized['Header']['time']*1e15))
            ax.set_xlabel(r'x [${\mu}\mathrm{m}$]')
            ax.legend()
        fig.suptitle(key)

        fig.tight_layout()
        fig.savefig(key.replace('/', '_') + '.png', dpi=320)

    def test_group_velocity(self):
        tx = {}
        for solver in self.solvers:
            tx[solver] = np.array([xt(dump)
                                   for dump in self.dumps[solver][1:]])
        print(tx)

        vg = dict(lehe_x=vg_lehe, lehe_custom=vg_lehe, optimized=vg_opt)

        for solver, data in tx.items():
            vg_sim = np.polyfit(data[:, 0], data[:, 1], 1)[0]

            # For reference, right here, right now the following line prints

            # yee 284766391.118 285957057.716 0.00416379510941
            # lehe_x 309981206.147 311627789.852 0.00528381536589
            # pukhov 291262060.412 292363351.796 0.00376685852363

            print('{:11} {:.3f} {:.3f} {:.12f}'.format(solver, vg_sim,
                  vg[solver], abs(vg_sim-vg[solver])/vg[solver]))

            assert np.isclose(vg_sim, vg[solver], rtol=0.01)


if __name__ == '__main__':
    unittest.main()
