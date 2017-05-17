#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the off module in antitile.
"""

import unittest
import doctest
from antitile import off

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(off))
    return tests

if __name__ == '__main__':
    unittest.main()
