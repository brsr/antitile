#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the flat module in antitile.
"""

import unittest
import doctest
from antitile import flat

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(flat))
    return tests

if __name__ == '__main__':
    unittest.main()
