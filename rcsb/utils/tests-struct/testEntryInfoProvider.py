##
# File:    EntryInfoProviderTests.py
# Author:  J. Westbrook
# Date:    22-Sep-2021
#
# Update:
#
#
##
"""
Tests for accessors (only) entry-level annotations.

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import unittest

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.struct.EntryInfoProvider import EntryInfoProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class EntryInfoProviderTests(unittest.TestCase):
    doInternal = True

    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__dataPath = os.path.join(HERE, "test-data")
        dirPath = os.path.join(self.__cachePath, "rcsb_entry_info")
        fU = FileUtil()
        fn = "entry_info_details.json"
        fU.put(os.path.join(self.__dataPath, fn), os.path.join(dirPath, fn))
        #

    def tearDown(self):
        pass

    def testGetEntryInfo(self):
        minCount = 12
        eiP = EntryInfoProvider(cachePath=self.__cachePath, useCache=True)
        ok = eiP.testCache(minCount=0)
        self.assertTrue(ok)
        riD = eiP.getEntryInfo("4en8")
        logger.info("riD (%d) %r", len(riD), riD)
        rL = eiP.getEntriesByPolymerEntityCount(count=2)
        self.assertGreaterEqual(len(rL), 5)
        ok = eiP.testCache(minCount=minCount)
        self.assertTrue(ok)


def entryInfoSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(EntryInfoProviderTests("testGetEntryInfo"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = entryInfoSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
