#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the xmath module in antitile.
"""

import unittest
import doctest
from antitile import xmath

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(xmath))
    return tests

if __name__ == '__main__':
    unittest.main()
