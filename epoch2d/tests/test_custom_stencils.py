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

import numpy as np
import sdf
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import os.path as osp
import unittest
import platform
from . import SimTest

micron = 1e-6
femto = 1e-15
c = 2.99792458e8 # m/s

# check that these correspond to input deck!

nx = 240
ny = nx//3

x_min = -12 * micron
x_max = -x_min
y_min = x_min
y_max = x_max

lambda_l = 0.5 * micron
x0 = -12.0 * micron # m
t0 = 8 * femto # s

# derived quantities from the above
k_l = 2*np.pi/lambda_l
dx = (x_max-x_min)/nx
kappax = k_l * dx
print('kappax: ', kappax)
dy = (y_max-y_min)/ny

#dt_lehe   = 0.95 * 1.0/np.sqrt(max(1.0/dx**2, 1.0/dy**2)) / c
dt_yee    = 0.99 * dx * dy / np.sqrt(dx**2+dy**2) / c



def showdatafields(sdffile):
    data = sdf.read(sdffile, dict=True)
    print('Data fields present in "' + sdffile + '":')
    print(data.keys())

def plotdump(sdffile, key, ax):
    data = sdf.read(sdffile, dict=True)
    array = data[key].data
    xaxis = data[key].grid_mid.data[0]*1e6
    yaxis = data[key].grid_mid.data[1]*1e6
    im = ax.imshow(array.T, extent=[min(xaxis), max(xaxis), min(yaxis), max(yaxis)])
    im.set_clim([-1e11,1e11])
    t = data['Header']['time']
    dx = np.asscalar(data[key].grid_mid.data[0][1] - data[key].grid_mid.data[0][0])

    #vg_lehe = c*(1.0 + 2.0*(1.0-c*dt_lehe/dx)*(k_l*dx/2.0)**2)
    #vg_yee  = c*(1.0 - (1.0-c*dt/dx)*(k_l*dx/2.0)**2)
    vg_yee = c*np.cos(k_l*dx/2.0)/np.sqrt(1-(c*dt_yee/dx*np.sin(k_l*dx/2.0))**2)

    # calculate_omega.py --dt-multiplier 0.99 --Y 3 --params 0.9082126568805592 0.04075757835916255
    # -0.04142032920970152 -0.20827814817872584 --physical-dx 1e-7 --laser-wavelength 5e-7
    # yields vg=1.0490493627815458 c
    vg_opt = 1.0490493627815458 * c

    # calculate_omega.py --dt-multiplier 0.99 --Y 3 --symmetric  --params 0.8988685682513151
    # 0.01862292597327679 -0.04155873453935287 --physical-dx 1e-7 --laser-wavelength 5e-7
    # yields vg=1.044753207834214 c
    vg_opt_symm = 1.044753207834214 * c

    # calculate_omega.py --dt-multiplier 0.99 --Y 3  --params 0.956632159129662  0.025096871992206993
    # -0.017744324957063393 -0.0009692545471922645 --physical-dx 1e-7 --laser-wavelength 5e-7
    # yields vg=1.0197513694119302 c
    vg_opt_xaxis = 1.0197513694119302 * c


    x = (x0 + c*(t-t0))*1e6
    x_opt = (x0 + vg_opt*(t-t0))*1e6
    x_opt_xaxis = (x0 + vg_opt_xaxis*(t-t0))*1e6
    x_yee = (x0 + vg_yee*(t-t0))*1e6

    if x > 0:
        ax.plot([x, x], [min(yaxis), 0.9*min(yaxis)], color='k', linestyle='-', linewidth=2, alpha=0.25)
        ax.plot([x, x], [0.9*max(yaxis), max(yaxis)], color='k', linestyle='-', linewidth=2, alpha=0.25)

        ax.plot([x_opt, x_opt], [min(yaxis), 0.9*min(yaxis)], color='b', linestyle='-', linewidth=2, alpha=0.25)
        ax.plot([x_opt, x_opt], [0.9*max(yaxis), max(yaxis)], color='b', linestyle='-', linewidth=2, alpha=0.25)

        ax.plot([x_opt_xaxis, x_opt_xaxis], [min(yaxis), 0.9*min(yaxis)], color='y', linestyle='-', linewidth=2, alpha=0.25)
        ax.plot([x_opt_xaxis, x_opt_xaxis], [0.9*max(yaxis), max(yaxis)], color='y', linestyle='-', linewidth=2, alpha=0.25)

        ax.plot([x_yee, x_yee], [min(yaxis), 0.9*min(yaxis)], color='r', linestyle='-', linewidth=2, alpha=0.25)
        ax.plot([x_yee, x_yee], [0.9*max(yaxis), max(yaxis)], color='r', linestyle='-', linewidth=2, alpha=0.25)

    ax.set_title('{:2.1f} fs'.format(data['Header']['time']*1e15))

def plotevolution(path, key):
    sdffiles = [osp.join(path, '{:04d}.sdf'.format(i)) for i in range(3)]

    fig, axarr = plt.subplots(1,3, figsize=(10,4))
    axarr[0].set_ylabel('y [µm]')
    for (sdffile, ax) in zip(sdffiles, np.ravel(axarr)):

        plotdump(sdffile, key, ax)
        ax.set_xlabel('x [µm]')

    if path.endswith('optimized'):
        title = 'Opt'
    elif path.endswith('optimized_symm'):
        title = 'Opt Symm'
    elif path.endswith('optimized_xaxis'):
        title = 'Opt xaxis'

    fig.suptitle(key)

    fig.tight_layout()
    fig.savefig(key.replace('/','_') + '_' + title + '.png', dpi=320)

def createplots(path):
    plotevolution(path, 'Electric Field/Ey')


class test_custom_stencils(SimTest):
    def __init__(self, arg):
        super(test_custom_stencils, self).__init__(arg)
        if platform.system() == 'Darwin':
            print('macosx backend')
            plt.switch_backend('macosx')

    def test_createplots(self):
        base = os.getcwd()

        createplots(osp.join(base, 'optimized'))
        createplots(osp.join(base, 'optimized_symm'))
        createplots(osp.join(base, 'optimized_xaxis'))



if __name__=='__main__':
    unittest.main()


