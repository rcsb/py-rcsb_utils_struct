##
#  File:           EntryInfoProvider.py
#  Date:           22-Sep-2021 jdw
#
#  Updated:
#
##
"""
Accessors for entry-level annotations extracted from ExDB.

"""

import logging
import os.path
import time

from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class EntryInfoProvider(StashableBase):
    """Accessors (only) for entry-level annotations."""

    def __init__(self, **kwargs):
        #
        self.__version = "0.50"
        cachePath = kwargs.get("cachePath", ".")
        useCache = kwargs.get("useCache", True)
        self.__dirName = "rcsb_entry_info"
        self.__dirPath = os.path.join(cachePath, self.__dirName)
        super(EntryInfoProvider, self).__init__(cachePath, [self.__dirName])
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__entryInfoD = self.__reload(fmt="json", useCache=useCache)
        #

    def testCache(self, minCount=1):
        if minCount == 0:
            return True
        if self.__entryInfoD and minCount and "entryInfo" in self.__entryInfoD and len(self.__entryInfoD["entryInfo"]) > minCount:
            logger.info("Entry annotations for (%d) entries", len(self.__entryInfoD["entryInfo"]))
            return True
        return False

    def getEntryInfo(self, entryId):
        """Return a dictionary of entry-level annotations.

        Returns:
            (dict): of entry-level annotations
        """
        try:
            return self.__entryInfoD["entryInfo"][entryId.upper()] if entryId.upper() in self.__entryInfoD["entryInfo"] else {}
        except Exception as e:
            logger.error("Failing with %r", str(e))
        return {}

    def getEntriesByPolymerEntityCount(self, count):
        oL = []
        try:
            for entryId, eD in self.__entryInfoD["entryInfo"].items():
                if eD["polymer_entity_count"] == count:
                    oL.append(entryId)
        except Exception as e:
            logger.error("Failing with %r", str(e))
        return oL

    def __getEntryInfoFilePath(self, fmt="json"):
        baseFileName = "entry_info_details"
        fExt = ".json" if fmt == "json" else ".pic"
        fp = os.path.join(self.__dirPath, baseFileName + fExt)
        return fp

    def reload(self):
        """Reload from the current cache file."""
        ok = False
        try:
            self.__entryInfoD = self.__reload(fmt="json", useCache=True)
            ok = self.__entryInfoD is not None
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __reload(self, fmt="json", useCache=True):
        entryInfoFilePath = self.__getEntryInfoFilePath(fmt=fmt)
        tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        pcD = {"version": self.__version, "created": tS, "identifiers": {}}

        if useCache and self.__mU.exists(entryInfoFilePath):
            logger.info("Reading entry-info cached path %r", entryInfoFilePath)
            pcD = self.__mU.doImport(entryInfoFilePath, fmt=fmt)
        return pcD
