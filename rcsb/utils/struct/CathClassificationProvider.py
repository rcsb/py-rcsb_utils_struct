##
#  File:  CathClassificationProvider.py
#  Date:  3-Apr-2019 jdw
#
#  Updates:
#  15-Jul-2021 jdw Update the constructor to common provider API conventions
##
"""
  Extract CATH domain assignments, term descriptions and CATH classification hierarchy
  from CATH flat files.

"""

import collections
import logging
import os.path
import sys
from datetime import datetime
from datetime import timedelta

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class CathClassificationProvider(StashableBase):
    """Extract CATH domain assignments, term descriptions and CATH classification hierarchy
    from CATH flat files.
    """

    def __init__(self, **kwargs):
        #
        self.__dirName = "cath"
        if "cachePath" in kwargs:
            self.__cachePath = os.path.abspath(kwargs.get("cachePath", None))
            self.__cathDirPath = os.path.join(self.__cachePath, self.__dirName)
        else:
            self.__cathDirPath = kwargs.get("cathDirPath", ".")
            self.__cachePath, self.__dirName = os.path.split(os.path.abspath(self.__cathDirPath))
        super(CathClassificationProvider, self).__init__(self.__cachePath, [self.__dirName])
        #
        useCache = kwargs.get("useCache", True)
        urlTarget = kwargs.get("cathTargetUrl", "http://download.cathdb.info/cath/releases/daily-release/newest")
        urlFallbackTarget = kwargs.get("cathTargetUrl", "http://download.cathdb.info/cath/releases/daily-release/archive")
        # no trailing /
        urlBackupPath = kwargs.get("cathUrlBackupPath", "https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets/master/fall_back/CATH")
        #
        self.__mU = MarshalUtil(workPath=self.__cathDirPath)
        self.__nD, self.__pdbD = self.__reload(urlTarget, urlFallbackTarget, self.__cathDirPath, useCache=useCache)
        if not self.testCache() and not useCache:
            ok = self.__fetchFromBackup(urlBackupPath, self.__cathDirPath)
            if ok:
                self.__nD, self.__pdbD = self.__reload(urlTarget, urlFallbackTarget, self.__cathDirPath, useCache=True)
        #

    def testCache(self):
        logger.info("CATH lengths nD %d pdbD %d", len(self.__nD), len(self.__pdbD))
        if (len(self.__nD) > 100) and (len(self.__pdbD) > 5000):
            return True
        return False

    def getCathVersions(self, pdbId, authAsymId):
        """aD[(pdbId, authAsymId)] = [(cathId, domainId, (authAsymId, resBeg, resEnd), version)]"""
        try:
            return list(set([tup[3] for tup in self.__pdbD[(pdbId, authAsymId)]]))
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))

        return []

    def getCathIds(self, pdbId, authAsymId):
        try:
            return list(set([tup[0] for tup in self.__pdbD[(pdbId, authAsymId)]]))
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))

        return []

    def getCathDomainNames(self, pdbId, authAsymId):
        try:
            return list(set([tup[1] for tup in self.__pdbD[(pdbId, authAsymId)]]))
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))

        return []

    def getCathResidueRanges(self, pdbId, authAsymId):
        try:
            return [(tup[0], tup[1], tup[2][0], tup[2][1], tup[2][2]) for tup in self.__pdbD[(pdbId, authAsymId)]]
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))

        return []

    def getCathName(self, cathId):
        try:
            return self.__nD[cathId]
        except Exception:
            logger.debug("Undefined CATH id %r", cathId)
        return None

    def getIdLineage(self, cathId):
        try:
            ff = cathId.split(".")
            return [".".join(ff[0:jj]) for jj in range(1, len(ff) + 1)]
        except Exception:
            logger.debug("No lineage for bad CATH id %r", cathId)
        return None

    def getNameLineage(self, cathId):
        try:
            return [self.getCathName(cId) for cId in self.getIdLineage(cathId)]
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return None

    def getTreeNodeList(self):
        return self.__exportTreeNodeList(self.__nD)

    def __getCathDomainFileName(self):
        pyVersion = sys.version_info[0]
        fn = "cath_domains-py%s.pic" % str(pyVersion)
        return fn

    def __reload(self, urlTarget, urlFallbackTarget, cathDirPath, useCache=True):
        nD = {}
        pdbD = {}
        fn = self.__getCathDomainFileName()
        cathDomainPath = os.path.join(cathDirPath, fn)
        self.__mU.mkdir(cathDirPath)
        #
        # cathDomainPath = os.path.join(cathDirPath, "cath_domains.json")
        #
        if useCache and self.__mU.exists(cathDomainPath):
            sD = self.__mU.doImport(cathDomainPath, fmt="pickle")
            logger.debug("Cath domain length %d", len(sD))
            nD = sD["names"]
            pdbD = sD["assignments"]
        elif not useCache:
            minLen = 1000
            logger.info("Fetch CATH name and domain assignment data from primary data source %s", urlTarget)
            nmL, dmL = self.__fetchFromSource(urlTarget, urlFallbackTarget, minLen)
            #
            ok = False
            nD = self.__extractNames(nmL)
            dD = self.__extractDomainAssignments(dmL)
            pdbD = self.__buildAssignments(dD)
            sD = {"names": nD, "assignments": pdbD}
            if (len(nD) > minLen) and (len(dD) > minLen):
                ok = self.__mU.doExport(cathDomainPath, sD, fmt="pickle")
            logger.debug("Cache save status %r", ok)
            #
        return nD, pdbD

    def __fetchFromBackup(self, urlBackupPath, cathDirPath):
        fn = self.__getCathDomainFileName()
        cathDomainPath = os.path.join(cathDirPath, fn)
        self.__mU.mkdir(cathDirPath)
        #
        backupUrl = urlBackupPath + "/" + fn
        logger.info("Using backup URL %r", backupUrl)
        fU = FileUtil()
        ok = fU.get(backupUrl, cathDomainPath)
        return ok

    def __fetchFromSource(self, urlTarget, urlFallbackTarget, minLen):
        """Fetch the classification names and domain assignments from CATH repo.

        http://download.cathdb.info/cath/releases/daily-release/newest/cath-b-newest-all.gz
        http://download.cathdb.info/cath/releases/daily-release/newest/cath-b-newest-names.gz
        #
        http://download.cathdb.info/cath/releases/daily-release/archive/cath-b-yyyymmdd-all.gz
        http://download.cathdb.info/cath/releases/daily-release/archive/cath-b-yyyymmdd-names-all.gz
        """
        fn = "cath-b-newest-names.gz"
        url = os.path.join(urlTarget, fn)
        nmL = self.__mU.doImport(url, fmt="list", uncomment=True)
        #
        if not nmL or len(nmL) < minLen:
            dS = datetime.today().strftime("%Y%m%d")
            dS = datetime.strftime(datetime.now() - timedelta(1), "%Y%m%d")
            fn = "cath-b-%s-names-all.gz" % dS
            url = os.path.join(urlFallbackTarget, fn)
            logger.info("Using fallback resource for %s", fn)
            nmL = self.__mU.doImport(url, fmt="list", uncomment=True)
        #
        fn = "cath-b-newest-all.gz"
        url = os.path.join(urlTarget, fn)
        dmL = self.__mU.doImport(url, fmt="list", uncomment=True)
        #
        if not dmL or len(dmL) < minLen:
            dS = datetime.today().strftime("%Y%m%d")
            dS = datetime.strftime(datetime.now() - timedelta(1), "%Y%m%d")
            fn = "cath-b-%s-all.gz" % dS
            url = os.path.join(urlFallbackTarget, fn)
            logger.info("Using fallback resource for %s", fn)
            dmL = self.__mU.doImport(url, fmt="list", uncomment=True)
        #
        return nmL, dmL

    def __extractNames(self, nmL):
        """
        From cath-b-newest-names:

            1 Mainly Alpha
            2 Mainly Beta
            3 Alpha Beta
            4 Few Secondary Structures
            1.10 Orthogonal Bundle
            1.20 Up-down Bundle
            1.25 Alpha Horseshoe
            1.40 Alpha solenoid
            1.50 Alpha/alpha barrel
            2.10 Ribbon
            2.20 Single Sheet
            2.30 Roll
            2.40 Beta Barrel
            2.50 Clam
            2.60 Sandwich
            2.70 Distorted Sandwich
            2.80 Trefoil
            2.90 Orthogonal Prism
            2.100 Aligned Prism
            2.102 3-layer Sandwich
        """
        rD = {}
        logger.info("length of input name list %d", len(nmL))
        for nm in nmL:
            ff = nm.split(" ")
            rD[ff[0]] = " ".join(ff[1:])
        return rD

    def __extractDomainAssignments(self, dmL):
        """
        From cath-b-newest-all:

            101mA00 v4_2_0 1.10.490.10 0-153:A
            102lA00 v4_2_0 1.10.530.40 1-162:A
            102mA00 v4_2_0 1.10.490.10 0-153:A
            103lA00 v4_2_0 1.10.530.40 1-162:A
            103mA00 v4_2_0 1.10.490.10 0-153:A
            104lA00 v4_2_0 1.10.530.40 1-162:A
            104lB00 v4_2_0 1.10.530.40 1-162:B
            104mA00 v4_2_0 1.10.490.10 1-153:A
            105mA00 v4_2_0 1.10.490.10 1-153:A
            106mA00 v4_2_0 1.10.490.10 0-153:A
            107lA00 v4_2_0 1.10.530.40 1-162:A
            107mA00 v4_2_0 1.10.490.10 0-153:A
            108lA00 v4_2_0 1.10.530.40 1-162:A
            108mA00 v4_2_0 1.10.490.10 0-153:A
            109lA00 v4_2_0 1.10.530.40 1-162:A
            109mA00 v4_2_0 1.10.490.10 0-153:A
            10gsA01 v4_2_0 3.40.30.10 2-78:A,187-208:A

            #
            Returns:

            dD[domainId] = (cathId, [(authAsymId, resBeg, resEnd), ...], version)

        """
        dD = {}
        logger.info("length of input domain assignment list %d", len(dmL))
        for dm in dmL:
            #
            try:
                ff = dm.split(" ")
                #
                rngL = ff[3].split(",")
                dmTupL = []
                for rng in rngL:
                    tL = rng.split(":")
                    rL = tL[0].split("-")
                    dmTupL.append((tL[1], rL[0], rL[1]))
                #
                dD[ff[0]] = (ff[2], dmTupL, ff[1])
            except Exception:
                logger.info("Failing for case %r: %r", ff, dm)
        return dD

    def __buildAssignments(self, dD):
        """
          Input internal data structure with domain assignments -

          dD[domainId] = (cathId, rangelist, version)

          Returns:

        =
           aD[(pdbId, authAsymId)] = [(cathId, domainId, (authAsymId, resBeg, resEnd), version)]
        """
        pdbD = {}
        for domId, dTup in dD.items():
            pdbId = domId[:4]
            for rTup in dTup[1]:
                pdbD.setdefault((pdbId, rTup[0]), []).append((dTup[0], domId, rTup, dTup[2]))
        return pdbD

    def __exportTreeNodeList(self, nD):
        """Create node list from name dictionary and lineage dictionaries."""
        # create parent dictionary
        #
        pL = []
        pD = {}
        for tId in nD:
            ff = tId.split(".")
            if len(ff) == 1:
                ptId = None
                pL.append(tId)
            else:
                ptId = ".".join(ff[:-1])
            logger.debug("tId %s parent %s", tId, ptId)
            pD[tId] = ptId
        #
        logger.info("nD %d pD %d", len(nD), len(pD))
        # create child dictionary
        cD = {}
        for ctId, ptId in pD.items():
            cD.setdefault(ptId, []).append(ctId)
        #
        logger.info("cD %d", len(cD))
        #
        idL = []
        for rootId in sorted(pL):
            visited = set([rootId])
            queue = collections.deque(visited)
            while queue:
                tId = queue.popleft()
                idL.append(tId)
                if tId not in cD:
                    # logger.debug("No children for CATH tId %s", tId)
                    continue
                for childId in cD[tId]:
                    if childId not in visited:
                        queue.append(childId)
                        visited.add(childId)
        #
        dL = []
        for tId in idL:
            displayName = nD[tId]
            ptId = pD[tId]
            ff = tId.split(".")
            lL = [".".join(ff[0:jj]) for jj in range(1, len(ff) + 1)]
            #
            # d = {'id': tId, 'name': displayName, 'lineage': lL, 'parents': [ptId], 'depth': len(lL)}
            if len(lL) == 1:
                dD = {"id": tId, "name": displayName, "depth": 0}
            else:
                dD = {"id": tId, "name": displayName, "parents": [ptId], "depth": len(lL) - 1}
            dL.append(dD)

        return dL
