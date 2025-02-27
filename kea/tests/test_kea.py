"""
Unit and regression test for the kea package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import kea


def test_kea_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "kea" in sys.modules
