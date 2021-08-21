##
#  File:  EcodClassificationProvider.py
#  Date:  3-Apr-2019 jdw
#
#  Updates:
#
##
"""
  Extract ECOD domain assignments, term descriptions and ECOD classification hierarchy
  from ECOD flat files.

"""

import collections
import datetime
import logging
import os.path
import sys

from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase

logger = logging.getLogger(__name__)


class EcodClassificationProvider(StashableBase):
    """Extract ECOD domain assignments, term descriptions and ECOD classification hierarchy
    from ECOD flat files.

    http://prodata.swmed.edu/ecod/

    See:
    H. Cheng, R. D. Schaeffer, Y. Liao, L. N. Kinch, J. Pei, S. Shi, B. H. Kim, N. V. Grishin. (2014)
    ECOD: An evolutionary classification of protein domains. PLoS Comput Biol 10(12): e1003926.

    Linking details:  http://prodata.swmed.edu/ecod/complete/domain/<domainId>

                      http://prodata.swmed.edu/ecod/complete/domain/e6sl5G1
    """

    #
    # --
    def __init__(self, cachePath, useCache, **kwargs):
        self.__cachePath = cachePath
        self.__useCache = useCache
        dirName = "ecod"
        super(EcodClassificationProvider, self).__init__(self.__cachePath, [dirName])
        self.__dirPath = os.path.join(cachePath, "ecod")
        self.__version = None
        #
        urlTarget = kwargs.get("ecodTargetUrl", "http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt")
        urlBackup = kwargs.get("ecodUrlBackupPath", "https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets/master/fall_back/ECOD/ecod.latest.domains.txt.gz")
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__pD, self.__nD, self.__ntD, self.__pdbD = self.__reload(urlTarget, urlBackup, self.__dirPath, useCache=useCache)

    def testCache(self):
        logger.info("ECOD Lengths nD %d pdbD %d", len(self.__nD), len(self.__pdbD))
        if (len(self.__nD) > 100) and (len(self.__pdbD) > 5000):
            return True
        return False

    def getVersion(self):
        return self.__version

    # --
    def getFamilyIds(self, pdbId, authAsymId):
        try:
            return list(set([tup[1] for tup in self.__pdbD[(pdbId.lower(), authAsymId)]]))
        except Exception as e:
            logger.exception("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getDomainIds(self, pdbId, authAsymId):
        try:
            return list(set([tup[0] for tup in self.__pdbD[(pdbId.lower(), authAsymId)]]))
        except Exception as e:
            logger.exception("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getFamilyNames(self, pdbId, authAsymId):
        try:
            return list(set([self.getName(tup[1]) for tup in self.__pdbD[(pdbId.lower(), authAsymId)]]))
        except Exception as e:
            logger.exception("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getFamilyResidueRanges(self, pdbId, authAsymId):
        try:
            # pdbD.setdefault((pdbId, authAsymId), []).append((domId, fId, authAsymId, authSeqBeg, authSeqEnd))
            return [(tup[0], tup[1], tup[2], tup[3], tup[4]) for tup in self.__pdbD[(pdbId.lower(), authAsymId)]]
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getName(self, domId):
        try:
            return self.__nD[domId].split("|")[0]
        except Exception:
            logger.debug("Undefined ECOD id %r", domId)
        return None

    def getNameType(self, domId):
        qD = {"A": "Architecture", "X": "Possible Homology", "H": "Homology", "T": "Topology", "F": "Family"}
        try:
            return qD[self.__ntD[domId]]
        except Exception:
            logger.debug("Undefined ECOD id %r", domId)
        return None

    def getIdLineage(self, domId):
        pList = []
        try:
            pList.append(domId)
            if domId == 0:
                return pList
            pt = self.__pD[domId]
            while (pt is not None) and (pt != 0):
                pList.append(pt)
                pt = self.__pD[pt]
        except Exception as e:
            logger.exception("Failing for %r with %s", domId, str(e))
        #
        pList.reverse()
        return pList

    def getNameLineage(self, domId):
        try:
            nL = []
            for dId in self.getIdLineage(domId):
                tN = self.getName(dId)
                tN = tN if tN else "Unnamed"
                nL.append(tN)
            return nL
        except Exception as e:
            logger.exception("Failing for %r with %s", domId, str(e))
        return None

    def getTreeNodeList(self):
        return self.__exportTreeNodeList(self.__pD)

    def __getDomainFileName(self):
        pyVersion = sys.version_info[0]
        fn = "ecod_domains-py%s.pic" % str(pyVersion)
        return fn

    def __reload(self, urlTarget, urlBackup, ecodDirPath, useCache=True):
        pD = nD = ntD = pdbD = {}
        fn = self.__getDomainFileName()
        ecodDomainPath = os.path.join(ecodDirPath, fn)
        self.__mU.mkdir(ecodDirPath)
        #
        if useCache and self.__mU.exists(ecodDomainPath):
            sD = self.__mU.doImport(ecodDomainPath, fmt="pickle")
            logger.debug("ECOD domain length %d", len(sD))
            nD = sD["names"]
            ntD = sD["nametypes"]
            pD = sD["parents"]
            pdbD = sD["assignments"]
            self.__version = sD["version"]
        elif not useCache:
            minLen = 1000
            logger.info("Fetch ECOD name and domain assignment data from primary data source %s", urlTarget)
            nmL = self.__fetchFromSource(urlTarget)
            if not nmL:
                nmL = self.__fetchFromSource(urlBackup)
            #
            logger.info("ECOD raw file length (%d)", len(nmL))
            ok = False
            pD, nD, ntD, pdbD = self.__extractDomainHierarchy(nmL)
            #
            tS = datetime.datetime.now().isoformat()
            vS = self.__version
            sD = {"version": vS, "created": tS, "names": nD, "nametypes": ntD, "parents": pD, "assignments": pdbD}
            if (len(nD) > minLen) and (len(pD) > minLen):
                ok = self.__mU.doExport(ecodDomainPath, sD, fmt="pickle")
            logger.debug("Cache save status %r", ok)
            #
        return pD, nD, ntD, pdbD

    def __fetchFromSource(self, urlTarget):
        """Fetch the classification names and domain assignments from the ECOD repo."""
        fU = FileUtil()
        fn = fU.getFileName(urlTarget)
        fp = os.path.join(self.__dirPath, fn)
        if not fU.exists(fp):
            fU.get(urlTarget, fp)
        #
        with open(fp, "r", encoding="utf-8") as ifh:
            line = ifh.readline()
            line = ifh.readline()
            line = ifh.readline()
            ff = line[:-1].split()
            self.__version = ff[-1]
        #
        nmL = self.__mU.doImport(fp, fmt="list", uncomment=True)
        fU.remove(fp)
        #
        return nmL

    def __extractDomainHierarchy(self, nmL):
        """
        #/data/ecod/database_versions/v280/ecod.develop280.domains.txt
        #ECOD version develop280
        #Domain list version 1.6
        #Grishin lab (http://prodata.swmed.edu/ecod)
        #uid	ecod_domain_id	manual_rep	f_id	pdb	chain	pdb_range	seqid_range	unp_acc	arch_name	x_name	h_name	t_name	f_name	asm_status	ligand
        002728551	e7d2xA1	AUTO_NONREP	1.1.1	7d2x	A	A:-3-183	A:20-206	NO_UNP	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	F_UNCLASSIFIED
        002728572	e7d5aA2	AUTO_NONREP	1.1.1	7d5a	A	A:-3-183	A:20-206	NO_UNP	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	F_UNCLASSIFIED
        002726563	e7b1eA1	AUTO_NONREP	1.1.1	7b1e	A	A:46P-183	A:14-199	NO_UNP	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	F_UNCLASSIFIED
        002726573	e7b1pA2	AUTO_NONREP	1.1.1	7b1p	A	A:47P-183	A:15-199	NO_UNP	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	F_UNCLASSIFIED
        """
        assignD = {}
        pD = {}
        ntD = {}
        hD = {}
        pIdD = {}
        nmD = {}
        #
        logger.info("Length of input ECOD name list %d", len(nmL))
        for nm in nmL:
            ff = nm.split("\t")
            # uId = ff[0]
            # ecodId is the linkable identifier -
            ecodId = ff[1]
            entryId = ff[4].lower()
            authAsymId = ff[5]
            resRange = ff[6]
            #
            #  There are no unique identifiers published for the internal elements of the hierarchy
            #   so these are assigned here similar to scop -   There are also many unnamed nodes
            #   that are conventionally filled in from the leaf levels of the tree...
            aGroupOrg = "A: " + ff[9].replace('"', "")
            xGroupOrg = "X: " + ff[10].replace('"', "")
            hGroupOrg = "H: " + ff[11].replace('"', "")
            tGroupOrg = "T: " + ff[12].replace('"', "")
            fGroupOrg = "F: " + ff[13].replace('"', "")
            if hGroupOrg == "H: NO_H_NAME":
                hGroupOrg = tGroupOrg + "|(NO_H)"
            if xGroupOrg == "X: NO_X_NAME":
                xGroupOrg = hGroupOrg + "|(NO_X)"
            #
            fGroupOrg = fGroupOrg if fGroupOrg != "F_UNCLASSIFIED" else "Unmapped domain of " + tGroupOrg
            #
            # Remove redundancy in names and assign unique ids
            #
            aGroup = aGroupOrg
            xGroup = xGroupOrg + "|" + aGroupOrg
            hGroup = hGroupOrg + "|" + xGroupOrg + "|" + aGroupOrg
            tGroup = tGroupOrg + "|" + hGroupOrg + "|" + xGroupOrg
            fGroup = fGroupOrg + "|" + tGroupOrg
            #
            hD.setdefault("A", set()).add(aGroup)
            hD.setdefault("X", set()).add(xGroup)
            hD.setdefault("H", set()).add(hGroup)
            hD.setdefault("T", set()).add(tGroup)
            hD.setdefault("F", set()).add(fGroup)
            aId = 100000 + len(hD["A"])
            xId = 200000 + len(hD["X"])
            hId = 300000 + len(hD["H"])
            tId = 400000 + len(hD["T"])
            fId = 500000 + len(hD["F"])
            #
            #
            if xGroup in pD and pD[xGroup] != aGroup:
                logger.error("skipping %r multiple parents for xGroup %r  %r and %r ", ecodId, xGroup, pD[xGroup], aGroup)
                continue
            #
            if hGroup in pD and pD[hGroup] != xGroup:
                logger.error("skipping %r multiple parents for hGroup %r  %r and %r ", ecodId, hGroup, pD[hGroup], xGroup)
                continue
            #
            if tGroup in pD and pD[tGroup] != hGroup:
                logger.error("skipping %r multiple parents for tGroup %r  %r and %r ", ecodId, tGroup, pD[tGroup], hGroup)
                continue
            #
            if fGroup in pD and pD[fGroup] != tGroup:
                logger.error("skipping %r multiple parents for fGroup %r  %r and %r ", ecodId, fGroup, pD[fGroup], tGroup)
                continue

            if xId in pIdD and pIdD[xId] != aId:
                logger.error("skipped %r multiple parents for xId %r  %r and %r ", ecodId, xId, pIdD[xId], aId)
            #
            if hId in pIdD and pIdD[hId] != xId:
                logger.error("skipped %r multiple parents for hId %r  %r and %r ", ecodId, hId, pIdD[hId], xId)
            #
            if tId in pIdD and pIdD[tId] != hId:
                logger.error("skipped %r multiple parents for tId %r  %r and %r ", ecodId, tId, pIdD[tId], hId)
            #
            if fId in pIdD and pIdD[fId] != tId:
                logger.error("skipped %r multiple parents for fId %r  %r and %r ", ecodId, fId, pIdD[fId], tId)

            #
            pIdD[aId] = 0
            pIdD[xId] = aId
            pIdD[hId] = xId
            pIdD[tId] = hId
            pIdD[fId] = tId
            #
            nmD[aId] = aGroupOrg
            nmD[xId] = xGroupOrg
            nmD[hId] = hGroupOrg
            nmD[tId] = tGroupOrg
            nmD[fId] = fGroupOrg
            #
            ntD[aId] = "A"
            ntD[xId] = "X"
            ntD[hId] = "H"
            ntD[tId] = "T"
            ntD[fId] = "F"
            rL = self.__parseRanges(resRange)
            assignD[(entryId, authAsymId)] = [(ecodId, fId, t[0], t[1], t[2]) for t in rL]
            #
        return pIdD, nmD, ntD, assignD

    def __parseRanges(self, rS):
        rL = []
        authAsymId = authSeqBeg = authSeqEnd = None
        try:
            tSL = rS.split(",")
            for tS in tSL:
                fL = tS.split(":")
                authAsymId = fL[0]
                rS = fL[1]
                if rS[0] == "-":
                    authSeqBeg = -int(rS[1:].split("-")[0])
                    authSeqEnd = int(rS[1:].split("-")[1])
                else:
                    authSeqBeg = int(rS.split("-")[0])
                    authSeqEnd = int(rS.split("-")[1])
            rL.append((authAsymId, authSeqBeg, authSeqEnd))
        except Exception:
            pass
        return rL

    def __exportTreeNodeList(self, pD):
        """Create node list from name dictionary and lineage dictionaries."""
        #
        rootId = 0
        pL = [rootId]
        #
        logger.info("pD %d pL %r", len(pD), pL)
        # --
        #
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
                    # logger.debug("No children for Ecod tId %s", tId)
                    continue
                for childId in cD[tId]:
                    if childId not in visited:
                        queue.append(childId)
                        visited.add(childId)
        #
        dL = []
        for tId in idL:
            displayName = self.getName(tId)
            ptId = pD[tId] if tId in pD else None
            lL = self.getIdLineage(tId)[1:]
            #
            if tId == rootId:
                continue
            elif ptId == rootId:
                dD = {"id": str(tId), "name": displayName, "depth": 0}
            else:
                dD = {"id": str(tId), "name": displayName, "parents": [str(ptId)], "depth": len(lL)}
            dL.append(dD)

        return dL
