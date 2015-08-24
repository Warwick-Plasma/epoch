#!/bin/env python

# Stephan Kuschel, 2015

import unittest

def main():
    runner = unittest.TextTestRunner()
    suite = unittest.TestSuite()
    import laser.test_laser
    suite.addTest(laser.test_laser.test_laser())
    print(suite)


if __name__=='__main__':
    main()
