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

def main():
    runner = unittest.TextTestRunner(verbosity=0)
    suite = getsuite('.')
    runner.run(suite)


if __name__=='__main__':
    main()
