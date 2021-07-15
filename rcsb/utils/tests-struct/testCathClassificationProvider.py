##
# File:    CathClassificationProviderTests.py
# Date:    3-Apr-2019  JDW
#
# Updates:
#
##
"""
Test cases for operations that read CATH term and class data from flat files -

"""

import logging
import os
import time
import unittest

from rcsb.utils.struct import __version__
from rcsb.utils.struct.CathClassificationProvider import CathClassificationProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class CathClassificationProviderTests(unittest.TestCase):
    def setUp(self):
        self.__workPath = os.path.join(HERE, "test-output")
        self.__cachePath = os.path.join(self.__workPath, "CACHE")
        #
        self.__startTime = time.time()
        logger.debug("Running tests on version %s", __version__)
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testGetCathData(self):
        """Load latest CATH data"""
        try:
            ccu = CathClassificationProvider(cachePath=self.__cachePath, useCache=False)
            pdbIdL = [("10gs", "A"), ("4hrt", "A"), ("4hrt", "C"), ("4hrt", "E"), ("4hrt", "G"), ("1jwn", "A"), ("1jwn", "B"), ("1jwn", "C"), ("1jwn", "D")]
            #
            for pdbTup in pdbIdL:
                cathids = ccu.getCathIds(pdbTup[0], pdbTup[1])
                domains = ccu.getCathDomainNames(pdbTup[0], pdbTup[1])
                ranges = ccu.getCathResidueRanges(pdbTup[0], pdbTup[1])
                versions = ccu.getCathVersions(pdbTup[0], pdbTup[1])
                logger.info("pdbId %r authAsymId %r cathids %r domains %r ranges %r versions %r", pdbTup[0], pdbTup[1], cathids, domains, ranges, versions)

            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testCathClassificationAccessMethods(self):
        """Test dictionary methods --"""
        try:
            ccu = CathClassificationProvider(cachePath=self.__cachePath, useCache=True)
            logger.info("Cath name for 1.10.490.10 %s ", ccu.getCathName("1.10.490.10"))
            nL = ccu.getTreeNodeList()
            logger.info("Node list length %d", len(nL))
            logger.info("Nodes %r", nL[:20])
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def readCathData():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(CathClassificationProviderTests("testGetCathData"))
    suiteSelect.addTest(CathClassificationProviderTests("testCathClassificationAccessMethods"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = readCathData()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
