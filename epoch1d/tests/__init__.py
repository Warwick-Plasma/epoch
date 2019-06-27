# Copyright (C) 2009-2019 University of Warwick
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

import unittest
import os
import subprocess


class SimTest(unittest.TestCase):
    '''
    The base class for all tests which need to run an EPOCH input.deck
    before starting with the tests.

    The Workin directory is deduced from the name of the test class:
    'test_laser' will use 'laser/' as working dir,
    'test_twostream' uses 'twostream/' and so on.
    '''

    @classmethod
    def setUpClass(cls):
        # remove 'test_' at the beginning of the class name
        simdir = cls.__name__[5:]
        # print('Changing to simdir: {}'.format(simdir))
        os.chdir(os.path.join(os.path.dirname(__file__), simdir))
        exitcode = subprocess.call('make', shell=True)
        cls.epochexitcode = exitcode
        if exitcode != 0:
            # that means the execution of 'make' returned an error
            os.chdir('..')

    def setUp(self):
        if self.epochexitcode:
            self.fail('running EPOCH errored (exitcode {:})'
                      .format(self.epochexitcode))

    @classmethod
    def tearDownClass(cls):
        if "NOCLEAN" not in os.environ:
            subprocess.call('make clean', shell=True)
        os.chdir('..')
