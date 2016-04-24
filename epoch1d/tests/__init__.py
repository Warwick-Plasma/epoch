# Stephan Kuschel, 2016

import unittest
import os, subprocess

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
        #print('Changing to simdir: {}'.format(simdir))
        os.chdir(os.path.join(os.path.dirname(__file__), simdir))
        exitcode = subprocess.call('make', shell=True)
        cls.epochexitcode = exitcode
        if exitcode != 0:
            # that means the execution of 'make' returned an error
            os.chdir('..')

    def setUp(self):
        if self.epochexitcode:
            self.fail('running EPOCH errored (exitcode {:})'.format(self.epochexitcode))

    @classmethod
    def tearDownClass(cls):
        os.chdir('..')
