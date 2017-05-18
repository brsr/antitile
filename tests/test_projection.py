#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests for the projection module in antitile.
"""

import unittest
import doctest
from antitile import projection

def load_tests(loader, tests, ignore):
    optionflags = (doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS |
                   doctest.IGNORE_EXCEPTION_DETAIL)    
    tests.addTests(doctest.DocTestSuite(projection, optionflags=optionflags))
    return tests

if __name__ == '__main__':
    unittest.main()
