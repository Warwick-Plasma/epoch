#!/bin/env python3

# Stephan Kuschel, 2015

import unittest
import os, subprocess

def getsuite(folder):
    files = os.listdir(folder)
    testfiles = [f.strip('.py') for f in files if f.startswith('test_')]
    loader = unittest.TestLoader()
    testcases = []
    for f in testfiles:
        mod = __import__(f)
        testcases.append(loader.loadTestsFromModule(mod))
    return unittest.TestSuite(testcases)

def run_tests():
    runner = unittest.TextTestRunner(verbosity=0)
    suite = getsuite('.')
    runner.run(suite)

def cleanup():
    for d in ['landau', 'laser', 'twostream']:
        os.chdir(d)
        subprocess.call(['make', 'clean'])
        os.chdir('..')

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--clean', '-c', action='store_true')
    args = parser.parse_args()

    if args.clean:
        cleanup()
    else:
        run_tests()



if __name__=='__main__':
    main()
