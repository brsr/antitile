#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the gcopoly module in antitile.
"""

import unittest
import doctest
from antitile import gcopoly

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(gcopoly))
    return tests

if __name__ == '__main__':
    unittest.main()
