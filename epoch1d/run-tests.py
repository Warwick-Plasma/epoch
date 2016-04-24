#!/usr/bin/env python

# Stephan Kuschel, 2015-2016

import nose
import subprocess
import os

_curfiledir = os.path.dirname(os.path.abspath(__file__))

def setcwd(relative=None):
    # reset current working dir to path of this file
    os.chdir(_curfiledir)
    if relative:
        os.chdir(relative)

def compileepoch():
    setcwd()
    exitcode = subprocess.call('make', shell=True)
    if exitcode != 0:
        print('compiling EPOCH errored (exitcode {})'.format(exitcode))
    return exitcode


def run_tests(args):
    noseargv = ['']
    if args.test:
        noseargv += ['tests.test_' + args.test]
    setcwd()
    testsok = nose.run(argv=noseargv)
    return testsok


def clean():
    '''
    clean the tests directory and all its subdirectoryies
    by calling 'make clean' in each of them
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
    parser = argparse.ArgumentParser()
    parser.add_argument('test', nargs='?', help='''
    run only a single test specified by its name, i.e. 'laser'
    ''')
    parser.add_argument('--clean', '-c', action='store_true', help=clean.__doc__)
    args = parser.parse_args()

    if args.clean:
        clean()
        exit()
    epochexitcode = compileepoch()
    if epochexitcode != 0:
        exit(epochexitcode)
    testsok = run_tests(args)
    exit(int(not testsok))


if __name__=='__main__':
    main()
