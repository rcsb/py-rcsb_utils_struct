##
# File:    EcodClassificationProviderTests.py
# Date:    23-Jun-2021  JDW
#
# Updates:
#
##
"""
Test cases for operations that read ECOD classification data from flat files -
"""

import logging
import os
import time
import unittest

from rcsb.utils.struct import __version__
from rcsb.utils.struct.EcodClassificationProvider import EcodClassificationProvider

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class EcodClassificationProviderTests(unittest.TestCase):
    def setUp(self):
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        #
        self.__startTime = time.time()
        logger.debug("Running tests on version %s", __version__)
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.debug("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testAGetEcodData(self):
        """Test ECOD data access methods"""
        try:
            ecodP = EcodClassificationProvider(cachePath=self.__cachePath, useCache=False)

            pdbIdL = [("4iqn", "C"), ("4gtl", "A"), ("4bt2", "A"), ("4hrt", "A"), ("1pca", "A"), ("1kn6", "A"), ("1jqg", "A")]
            #
            ver = ecodP.getVersion()
            self.assertTrue(ver is not None)
            logger.info("ECOD version %r", ver)
            #
            # e4b2tQ1	AUTO_NONREP	2487.1.1.35	4b2t	Q	Q:219-248,Q:276-369	Q:219-248,Q:276-36
            for pdbTup in pdbIdL:
                fIds = ecodP.getFamilyIds(pdbTup[0], pdbTup[1])
                for fId in fIds:
                    nm = ecodP.getName(fId)
                    self.assertTrue(nm)
                    nmT = ecodP.getNameType(fId)
                    self.assertEqual(nmT, "Family")
                    logger.debug("pdbTup %r %r fId: %r name: %r type: %r", pdbTup[0], pdbTup[1], fId, nmT, nm)
                    #
                    lL = ecodP.getIdLineage(fId)
                    self.assertTrue(lL)
                    nL = ecodP.getNameLineage(fId)
                    self.assertTrue(nL)
                    logger.debug("pdbTup %r %r lineage ids   %r", pdbTup[0], pdbTup[1], lL)
                    logger.debug("pdbTup %r %r lineage names %r", pdbTup[0], pdbTup[1], nL)
                fns = ecodP.getFamilyNames(pdbTup[0], pdbTup[1])
                fRanges = ecodP.getFamilyResidueRanges(pdbTup[0], pdbTup[1])
                logger.info("pdbTup %r %r - %r %r %r", pdbTup[0], pdbTup[1], fIds, fns, fRanges)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def testEcodTreeMethods(self):
        """Test ECOD tree methods --"""
        try:
            ccu = EcodClassificationProvider(cachePath=self.__cachePath, useCache=True)
            nL = ccu.getTreeNodeList()
            logger.info("ECOD node list length %d", len(nL))
            self.assertGreaterEqual(len(nL), 20000)
            logger.info("ECOD nodes %r", nL[:60])
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def ecodProviderSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(EcodClassificationProviderTests("testAGetEcodData"))
    suiteSelect.addTest(EcodClassificationProviderTests("testEcodTreeMethods"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = ecodProviderSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
