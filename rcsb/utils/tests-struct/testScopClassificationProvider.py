##
# File:    ScopClassificationProviderTests.py
# Date:    3-Apr-2019  JDW
#
# Updates:
#
##
"""
Test cases for operations that read SCOP term and class data from flat files -

"""

import logging
import os
import time
import unittest

from rcsb.utils.struct import __version__
from rcsb.utils.struct.ScopClassificationProvider import ScopClassificationProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class ScopClassificationProviderTests(unittest.TestCase):
    def setUp(self):
        self.__dirPath = os.path.join(os.path.dirname(TOPDIR), "rcsb", "mock-data")
        self.__workPath = os.path.join(HERE, "test-output")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__startTime = time.time()
        logger.debug("Running tests on version %s", __version__)
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testGetScopData(self):
        """Load latest scop data and test accessors"""
        try:
            scu = ScopClassificationProvider(cachePath=self.__cachePath, useCache=False)
            self.assertEqual(scu.getScopName(58788), "Designed proteins")
            self.assertEqual(scu.getNameLineage(58231), ["Peptides"])
            self.assertEqual(scu.getIdLineage(58231), [58231])
            #
            sunIdL = [46456, 48724, 51349, 53931, 56572, 56835, 56992, 57942, 58117, 58231, 58788, 310555]
            for sunId in sunIdL:
                logger.debug("ScopId %r name %s", sunId, scu.getScopName(sunId))
                logger.debug("ScopId %r lineage %r", sunId, scu.getNameLineage(sunId))
                logger.debug("ScopId %r lineage %r", sunId, scu.getIdLineage(sunId))
            #
            pdbIdL = [("4hrt", "A"), ("4hrt", "C"), ("4hrt", "E"), ("4hrt", "G"), ("1jwn", "A"), ("1jwn", "B"), ("1jwn", "C"), ("1jwn", "D")]
            #
            for pdbTup in pdbIdL:
                sunids = scu.getScopSunIds(pdbTup[0], pdbTup[1])
                domains = scu.getScopDomainNames(pdbTup[0], pdbTup[1])
                ranges = scu.getScopResidueRanges(pdbTup[0], pdbTup[1])
                logger.debug("pdbId %r authAsymId %r sunids %r domains %r ranges %r", pdbTup[0], pdbTup[1], sunids, domains, ranges)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testScopClassificationAccessMethods(self):
        """Test Scop tree node list generation --"""
        try:
            scu = ScopClassificationProvider(cachePath=self.__cachePath, useCache=True)
            nL = scu.getTreeNodeList()
            logger.debug("Node list length %d", len(nL))
            logger.debug("Nodes %r", nL[:30])
            self.assertGreaterEqual(len(nL), 22100)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def readScopData():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ScopClassificationProviderTests("testGetScopData"))
    suiteSelect.addTest(ScopClassificationProviderTests("testScopClassificationAccessMethods"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = readScopData()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
