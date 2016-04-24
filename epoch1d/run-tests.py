#!/bin/env python3

# Stephan Kuschel, 2015-2016

import nose
import subprocess
import os

def setcwd():
    os.chdir(os.path.dirname(__file__))

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


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('test', nargs='?')
    args = parser.parse_args()

    epochexitcode = compileepoch()
    if epochexitcode != 0:
        exit(epochexitcode)
    testsok = run_tests(args)
    exit(int(not testsok))


if __name__=='__main__':
    main()
