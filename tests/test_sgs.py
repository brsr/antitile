#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the sgs module in antitile.
"""

import unittest
import doctest
from antitile import sgs

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(sgs))
    return tests

if __name__ == '__main__':
    unittest.main()
