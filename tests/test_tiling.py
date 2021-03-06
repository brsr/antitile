#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the tiling module in antitile.
"""

import unittest
import doctest
from antitile import tiling

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(tiling))
    return tests

if __name__ == '__main__':
    unittest.main()
