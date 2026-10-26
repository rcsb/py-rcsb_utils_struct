##
# File:    CathClassificationProviderTests.py
# Date:    3-Apr-2019  JDW
#
# Updates:
#  26-Oct-2026 dwp Add test for getting entry IDs from short or extended PDB ID domain IDs, and fail testGetCathData
#                  once the CATH source data uses extended PDB IDs
#
##
"""
Test cases for operations that read CATH term and class data from flat files -

"""

import logging
import os
import time
import unittest

from importlib.metadata import version as get_package_version
from rcsb.utils.struct.CathClassificationProvider import CathClassificationProvider

__version__ = get_package_version("rcsb.utils.struct")

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(HERE))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


def getNonShortEntryIds(entryIds):
    """Return the entry IDs that aren't short (4-character) PDB IDs, e.g. extended (12-character) PDB IDs."""
    return sorted({eId for eId in entryIds if len(eId) != 4})


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

    def testEntryIdFromDomainId(self):
        """Get the lower case entry ID from short or extended PDB ID domain IDs (no data download)"""
        getEntryId = CathClassificationProvider._CathClassificationProvider__getEntryIdFromDomainId  # pylint: disable=protected-access
        for domId, entryId in [("101mA00", "101m"), ("10GSA01", "10gs"), ("1abc000", "1abc"), ("pdb_0000101mA00", "pdb_0000101m"), ("PDB_00001ABCB02", "pdb_00001abc")]:
            self.assertEqual(getEntryId(domId), entryId)

    # @unittest.skipIf(True, "Service temporarily down")
    def testGetCathData(self):
        """Load latest CATH data"""
        try:
            ccu = CathClassificationProvider(cachePath=self.__cachePath, useCache=False)
            pdbIdL = [
                ("10gs", "A"),
                ("4hrt", "A"),
                ("4hrt", "C"),
                ("4hrt", "E"),
                ("4hrt", "G"),
                ("1jwn", "A"),
                ("1jwn", "B"),
                ("1jwn", "C"),
                ("1jwn", "D"),
                ('3ehu', 'B'),
                ('3n93', 'B'),
                ('4qsk', 'B'),
            ]
            #
            for pdbTup in pdbIdL:
                cathids = ccu.getCathIds(pdbTup[0], pdbTup[1])
                domains = ccu.getCathDomainNames(pdbTup[0], pdbTup[1])
                ranges = ccu.getCathResidueRanges(pdbTup[0], pdbTup[1])
                versions = ccu.getCathVersions(pdbTup[0], pdbTup[1])
                logger.info("pdbId %r authAsymId %r cathids %r domains %r ranges %r versions %r", pdbTup[0], pdbTup[1], cathids, domains, ranges, versions)
                # Check for duplicate data items
                self.assertTrue(len(ranges) == len(set(ranges)))
            #
            # Fails once CATH switches to extended PDB IDs (also checks the domain ID length, <entryId><chain code><2-digit domain number>)
            pdbD = ccu._CathClassificationProvider__pdbD  # pylint: disable=protected-access
            badL = getNonShortEntryIds({pdbId for pdbId, _ in pdbD})
            badL += sorted({tup[1] for tupL in pdbD.values() for tup in tupL if len(tup[1]) != 7})
            self.assertEqual(badL, [], "CATH source data has entry IDs that aren't 4-char PDB IDs (update struct and dictionary for extended IDs): %r" % badL[:10])
        except AssertionError:
            raise

            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    # @unittest.skipIf(True, "Service temporarily down")
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
    suiteSelect.addTest(CathClassificationProviderTests("testEntryIdFromDomainId"))
    suiteSelect.addTest(CathClassificationProviderTests("testGetCathData"))
    suiteSelect.addTest(CathClassificationProviderTests("testCathClassificationAccessMethods"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = readCathData()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
