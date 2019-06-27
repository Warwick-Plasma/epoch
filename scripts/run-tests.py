#!/usr/bin/env python

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

import nose
import subprocess
import os

_epochdir = os.path.dirname(os.path.abspath(__file__))
_epochdir = os.path.join(_epochdir, '..')
_subdir = None


def setcwd(relative=None):
    '''
    resets the current working directiory to the path
    of this file.
    '''
    os.chdir(_epochdir)
    os.chdir(_subdir)
    if relative:
        os.chdir(relative)


def compileepoch():
    '''
    compiles the EPOCH code. The exit code of 'make' is returned.
    '''
    setcwd()
    exitcode = subprocess.call('make', shell=True)
    if exitcode != 0:
        print('compiling EPOCH errored (exitcode {})'.format(exitcode))
    return exitcode


def run_tests(args):
    '''
    use nose to collect the tests and run them all.
    '''
    noseargv = ['']
    if args.test:
        noseargv += ['tests.test_' + args.test]
    setcwd()
    testsok = nose.run(argv=noseargv)
    return testsok


def clean():
    '''
    clean the tests directory and all its subdirectoryies
    by calling 'make clean' in each of them.
    '''
    setcwd()
    subprocess.call('rm -rf tests/__pycache__', shell=True)  # python3
    subprocess.call('rm -rf tests/*.pyc', shell=True)  # python2
    setcwd()
    files = [os.path.join('tests', f) for f in os.listdir('tests')]
    dirs = [d for d in files if os.path.isdir(d)]
    for d in dirs:
        # call 'make clean' in every subdir
        setcwd(d)
        subprocess.call('make clean', shell=True)


def main():
    import argparse
    parser = argparse.ArgumentParser(description='''
    This runs the testsuite for EPOCH1D.
    It compiles EPOCH and runs the tests.
    It does NOT: install the python SDF reader or any other dependencies,
    which might be needed!
    ''')
    parser.add_argument('codeversion', help='''
    specify the  code version to run the tests on.
    ''', choices=['epoch1d', '1d', 'epoch2d', '2d', 'epoch3d', '3d'])
    parser.add_argument('test', nargs='?', help='''
    run only a single test specified by its name, i.e. 'laser'
    ''')
    parser.add_argument('--clean', '-c', action='store_true',
                        help=clean.__doc__)
    parser.add_argument('--build', '-b', action='store_true', help='''
    build only. Do not run the code.
    ''')
    args = parser.parse_args()
    subdirdict = {
        '1d': 'epoch1d',
        '2d': 'epoch2d',
        '3d': 'epoch3d'}
    global _subdir
    _subdir = args.codeversion
    if args.codeversion in subdirdict:
        _subdir = subdirdict[_subdir]

    if args.clean:
        clean()
        exit()
    epochexitcode = compileepoch()
    if epochexitcode != 0:
        exit(epochexitcode)
    if args.build:
        exit()
    testsok = run_tests(args)
    exit(int(not testsok))


if __name__ == '__main__':
    main()
