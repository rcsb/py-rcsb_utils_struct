##
# File:    Scop2ClassificationProviderTests.py
# Date:    23-Jun-2021  JDW
#
# Updates:
#
##
"""
Test cases for operations that read SCOP2 term and class data from flat files -
"""

import logging
import os
import time
import unittest

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.struct import __version__
from rcsb.utils.struct.Scop2ClassificationProvider import Scop2ClassificationProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class Scop2ClassificationProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__treePath = os.path.join(self.__cachePath, "scop2-tree.json")
        #
        self.__startTime = time.time()
        logger.debug("Running tests on version %s", __version__)
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testAGetScop2Data(self):
        """Test SCOP2 data access methods"""
        try:
            scp = Scop2ClassificationProvider(cachePath=self.__cachePath, useCache=False)

            pdbIdL = [("1pyt", "A"), ("4hrt", "A"), ("1pca", "A"), ("1kn6", "A"), ("1jqg", "A")]
            #
            ver = scp.getVersion()
            self.assertTrue(ver is not None)
            logger.info("SCOP2 version %r", ver)
            #
            for pdbTup in pdbIdL:
                fIds = scp.getFamilyIds(pdbTup[0], pdbTup[1])
                self.assertTrue(fIds)
                for fId in fIds:
                    nm = scp.getName(fId)
                    self.assertTrue(nm)
                    lL = scp.getIdLineage(fId)
                    self.assertTrue(lL)
                    nL = scp.getNameLineage(fId)
                    self.assertTrue(nL)
                    logger.debug("pdbTup %r %r lineage ids   %r", pdbTup[0], pdbTup[1], lL)
                    logger.debug("pdbTup %r %r lineage names %r", pdbTup[0], pdbTup[1], nL)
                sfIds = scp.getSuperFamilyIds(pdbTup[0], pdbTup[1])
                for sfId in sfIds:
                    lL = scp.getIdLineage(sfId)
                    self.assertTrue(lL)
                    nL = scp.getNameLineage(sfId)
                    self.assertTrue(nL)
                    logger.debug("pdbTup %r %r lineage ids   %r", pdbTup[0], pdbTup[1], lL)
                    logger.debug("pdbTup %r %r lineage names %r", pdbTup[0], pdbTup[1], nL)
                #
                fns = scp.getFamilyNames(pdbTup[0], pdbTup[1])
                self.assertTrue(fns)
                sfns = scp.getSuperFamilyNames(pdbTup[0], pdbTup[1])
                self.assertTrue(sfns)
                #
                fRanges = scp.getFamilyResidueRanges(pdbTup[0], pdbTup[1])
                self.assertTrue(fRanges)
                sfRanges = scp.getSuperFamilyResidueRanges(pdbTup[0], pdbTup[1])
                self.assertTrue(sfRanges)
                logger.debug("pdbTup %r %r - %r %r %r %r %r %r", pdbTup[0], pdbTup[1], fIds, sfIds, fns, sfns, fRanges, sfRanges)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testGetScop2BData(self):
        """Test SCOP2B data access methods"""
        try:
            scp = Scop2ClassificationProvider(cachePath=self.__cachePath, useCache=True)

            pdbIdL = [("5oba", "B"), ("6s61", "X")]
            #
            ver = scp.getVersion()
            self.assertTrue(ver is not None)
            logger.info("SCOP2 version %r", ver)
            #
            for pdbTup in pdbIdL:
                sfIds = scp.getSuperFamilyIds2B(pdbTup[0], pdbTup[1])
                self.assertTrue(sfIds)
                for sfId in sfIds:
                    lL = scp.getIdLineage(sfId)
                    self.assertTrue(lL)
                    nL = scp.getNameLineage(sfId)
                    self.assertTrue(nL)
                    logger.info("pdbTup %r %r lineage ids   %r", pdbTup[0], pdbTup[1], lL)
                    logger.info("pdbTup %r %r lineage names %r", pdbTup[0], pdbTup[1], nL)
                #
                sfns = scp.getSuperFamilyNames2B(pdbTup[0], pdbTup[1])
                sfRanges = scp.getSuperFamilyResidueRanges2B(pdbTup[0], pdbTup[1])
                logger.debug("pdbTup %r %r - %r %r %r", pdbTup[0], pdbTup[1], sfIds, sfns, sfRanges)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testScop2TreeMethods(self):
        """Test SCOP2 tree methods --"""
        try:
            ccu = Scop2ClassificationProvider(cachePath=self.__cachePath, useCache=True)
            nL = ccu.getTreeNodeList()
            self.assertGreaterEqual(len(nL), 75000)
            logger.info("SCOP2 tree node list length %d", len(nL))
            logger.debug("SCOP2 tree Nodes %r", nL[:20])
            mU = MarshalUtil(workPath=self.__cachePath)
            mU.doExport(self.__treePath, nL, fmt="json", indent=3)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def scop2ProviderSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(Scop2ClassificationProviderTests("testAGetScop2Data"))
    suiteSelect.addTest(Scop2ClassificationProviderTests("testGetScop2BData"))
    suiteSelect.addTest(Scop2ClassificationProviderTests("testScop2TreeMethods"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = scop2ProviderSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
