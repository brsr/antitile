#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the breakdown module in antitile.
"""

import unittest
import doctest
from antitile import breakdown

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(breakdown))
    return tests

if __name__ == '__main__':
    unittest.main()
