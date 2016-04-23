#!/bin/env python3

# Stephan Kuschel, 2015-2016

import nose
import subprocess
import os

def setcwd():
    os.chdir(os.path.dirname(__file__))

def compileepoch(exitonerror=True):
    setcwd()
    exitcode = subprocess.call('make', shell=True)
    if exitcode != 0:
        print('compiling EPOCH errored (exitcode {})'.format(exitcode))
        if exitonerror:
            exit(exitcode)
    return exitcode


def run_tests(args):
    compileepoch()
    noseargv = ['']
    if args.test:
        noseargv += ['tests.test_' + args.test]
    setcwd()
    nose.run(argv=noseargv)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('test', nargs='?')
    args = parser.parse_args()

    run_tests(args)



if __name__=='__main__':
    main()
