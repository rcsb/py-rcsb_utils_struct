##
#  File:  Scop2ClassificationProvider.py
#  Date:  29-Jun-2021 jdw
#
#  Updates:
#   10-Sep-2021 jdw split tree with type and class roots
#   24-Feb-2022 dwp Resolve duplication issues with Scop2 tree node list, and fix parent ID lists for nodes with multiple parents
#   18-Apr-2023 aae Use "pickle" as default file format
#   18-Jul-2023 dwp Resolve duplication issues with Scop2 families list
#   23-Apr-2024 dwp SCOP2/SCOP2B website was shut down--turn off fetching of source data until/if new site is made available again
##
"""
  Extract SCOP2 domain assignments, term descriptions and SCOP2 classification hierarchy
  from SCOP2 and SCOP2B flat files.

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


class Scop2ClassificationProvider(StashableBase):
    """Extract SCOP2 domain assignments, term descriptions and SCOP classification hierarchy
    from SCOP and SCOP2B flat files.
    """

    def __init__(self, cachePath, useCache, **kwargs):
        #
        _ = kwargs
        self.__cachePath = cachePath
        dirName = "scop2"
        self.__dirPath = os.path.join(self.__cachePath, dirName)
        self.__useCache = useCache
        super(Scop2ClassificationProvider, self).__init__(self.__cachePath, [dirName])
        #
        self.__version = "latest"
        self.__fmt = "pickle"
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        #
        # Temporarily turn off fetching of source data until new site is made available again
        # self.__nD, self.__ntD, self.__pAD, self.__pBD, self.__pBRootD, self.__fD, self.__sfD, self.__sf2bD = self.__reload(useCache=self.__useCache, fmt=self.__fmt)
        # if not useCache and not self.testCache():
        self.__fetchFromBackup(fmt=self.__fmt)
        self.__nD, self.__ntD, self.__pAD, self.__pBD, self.__pBRootD, self.__fD, self.__sfD, self.__sf2bD = self.__reload(useCache=True, fmt=self.__fmt)
        if not self.testCache():
            logger.error("Failed to build SCOP2 CACHE")

    def testCache(self):
        logger.info(
            "SCOP2 lengths nD %d pAD %d pBD %d pBRootD %d fD %d sfD %d sf2bD %d",
            len(self.__nD), len(self.__pAD), len(self.__pBD), len(self.__pBRootD), len(self.__fD), len(self.__sfD), len(self.__sf2bD)
        )
        if (len(self.__nD) > 9000) and (len(self.__pAD) > 70000):
            return True
        return False

    def getVersion(self):
        """Returns the SCOP2 version"""
        return self.__version

    def getFamilyIds(self, pdbId, authAsymId):
        try:
            return list(set([tup[1] for tup in self.__fD[(pdbId.upper(), authAsymId)]]))
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getSuperFamilyIds(self, pdbId, authAsymId):
        try:
            return list(set([tup[1] for tup in self.__sfD[(pdbId.upper(), authAsymId)]]))
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getFamilyNames(self, pdbId, authAsymId):
        try:
            return list(set([self.__nD[tup[1]] for tup in self.__fD[(pdbId.upper(), authAsymId)]]))
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getSuperFamilyNames(self, pdbId, authAsymId):
        try:
            return list(set([self.__nD[tup[1]] for tup in self.__sfD[(pdbId.upper(), authAsymId)]]))
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getFamilyResidueRanges(self, pdbId, authAsymId):
        try:
            # s/fD.setdefault((pdbId, authAsymId), []).append((domSuperFamilyId, authAsymId, authSeqBeg, authSeqEnd))
            return [(tup[0], tup[1], tup[2], tup[3], tup[4]) for tup in self.__fD[(pdbId.upper(), authAsymId)]]
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getSuperFamilyResidueRanges(self, pdbId, authAsymId):
        try:
            return [(tup[0], tup[1], tup[2], tup[3], tup[4]) for tup in self.__sfD[(pdbId.upper(), authAsymId)]]
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getSuperFamilyNames2B(self, pdbId, authAsymId):
        try:
            return list(set([self.__nD[tup[1]] for tup in self.__sf2bD[(pdbId.upper(), authAsymId)]]))
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getSuperFamilyIds2B(self, pdbId, authAsymId):
        try:
            return list(set([tup[1] for tup in self.__sf2bD[(pdbId.upper(), authAsymId)]]))
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getSuperFamilyResidueRanges2B(self, pdbId, authAsymId):
        try:
            return [(tup[0], tup[1], tup[2], tup[3], tup[4]) for tup in self.__sf2bD[(pdbId.upper(), authAsymId)]]
        except Exception as e:
            logger.debug("Failing for %r %r with %s", pdbId, authAsymId, str(e))
        return []

    def getName(self, domId):
        try:
            return self.__nD[domId]
        except Exception:
            logger.debug("Undefined SCOP2 id %r", domId)
        return None

    def getNameType(self, domId):
        qD = {"TP": "Protein Type", "CL": "Protein Class", "CF": "Fold", "SF": "Superfamily", "FA": "Family"}
        try:
            return qD[self.__ntD[domId]]
        except Exception:
            logger.debug("Undefined ECOD id %r", domId)
        return None

    def getIdLineage(self, domId):
        pS = set()
        try:
            pS.add(domId)
            pt = self.__pAD[domId]
            while (pt is not None) and (pt != 0):
                pS.add(pt)
                pt = self.__pAD[pt]
            #
            pt = self.__pBD[domId]
            while (pt is not None) and (pt != 0):
                pS.add(pt)
                pt = self.__pBD[pt]
        except Exception as e:
            logger.debug("Failing for %r with %s", domId, str(e))
        #
        return sorted(pS)

    def getNameLineage(self, domId):
        try:
            nL = []
            for dId in self.getIdLineage(domId):
                tN = self.getName(dId)
                tN = tN if tN else "Unnamed"
                nL.append(tN)
            return nL
        except Exception as e:
            logger.debug("Failing for %r with %s", domId, str(e))
        return None

    def getTreeNodeList(self):
        tnL = self.__exportTreeNodeList(self.__nD, self.__pAD, self.__pBRootD)
        return tnL

    def __getAssignmentFileName(self, fmt="pickle"):
        ext = "json" if fmt == "json" else "pic"
        fn = "scop2_domain_assignments.%s" % ext
        return fn

    def __reload(self, useCache=True, fmt="pickle"):
        nD = ntD = pAD = pBD = pBRootD = fD = sfD = sf2bD = {}
        fn = self.__getAssignmentFileName(fmt=fmt)
        assignmentPath = os.path.join(self.__dirPath, fn)
        self.__mU.mkdir(self.__dirPath)
        #
        if useCache and self.__mU.exists(assignmentPath):
            sD = self.__mU.doImport(assignmentPath, fmt=fmt)
            logger.debug("Domain name count %d", len(sD["names"]))
            self.__version = sD["version"]
            nD = sD["names"]
            ntD = sD["nametypes"]
            pAD = sD["parentsType"]
            pBD = sD["parentsClass"]
            pBRootD = sD["parentsClassRoot"]
            fD = sD["families"]
            sfD = sD["superfamilies"]
            sf2bD = sD["superfamilies2b"]

        elif not useCache:
            nmL, dmL, scop2bL, _ = self.__fetchFromSource()
            #
            ok = False
            nD = self.__extractNames(nmL)
            logger.info("Domain name dictionary (%d)", len(nD))
            pAD, pBD, pBRootD, ntD, fD, sfD, domToSfD = self.__extractDomainHierarchy(dmL)
            #
            logger.info("Domain node parent hierarchy (protein type) (%d)", len(pAD))
            logger.info("Domain node parent hierarchy (structural class) (%d)", len(pBD))
            logger.info("Domain node parent hierarchy (structural class root) (%d)", len(pBRootD))
            logger.info("SCOP2 core domain assignments (family %d) (sf %d)", len(fD), len(sfD))
            #
            sf2bD = self.__extractScop2bSuperFamilyAssignments(scop2bL, domToSfD)
            logger.info("SCOP2B SF domain assignments (%d)", len(sf2bD))
            #
            tS = datetime.datetime.now().isoformat()
            # vS = datetime.datetime.now().strftime("%Y-%m-%d")
            vS = self.__version
            sD = {
                "version": vS,
                "created": tS,
                "names": nD,
                "nametypes": ntD,
                "parentsType": pAD,
                "parentsClass": pBD,
                "parentsClassRoot": pBRootD,
                "families": fD,
                "superfamilies": sfD,
                "superfamilies2b": sf2bD
            }
            ok = self.__mU.doExport(assignmentPath, sD, fmt=fmt, indent=3)
            logger.info("Cache save status %r", ok)
            #
        return nD, ntD, pAD, pBD, pBRootD, fD, sfD, sf2bD

    def __fetchFromBackup(self, fmt="pickle"):
        urlTarget = "https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets/master/fall_back/SCOP2"
        #
        fn = self.__getAssignmentFileName(fmt=fmt)
        assignmentPath = os.path.join(self.__dirPath, fn)
        urlPath = os.path.join(urlTarget, fn)
        #
        logger.info("Creating directory %r", self.__dirPath)
        self.__mU.mkdir(self.__dirPath)
        logger.info("Fetching backup URL %r to local path %r", urlPath, assignmentPath)
        fU = FileUtil()
        ok = fU.get(urlPath, assignmentPath)
        logger.info("Fetch status %r", ok)
        return ok

    def __fetchFromSource(self):
        """Fetch the classification names and domain assignments from SCOP2 and SCOP2B resources.

        SCOP2 domain names:
            https://scop.mrc-lmb.cam.ac.uk/files/scop-des-latest.txt

        SCOP2 domain hierarchy:
            https://scop.mrc-lmb.cam.ac.uk/files/scop-cla-latest.txt

        SIFTS extrapolated SCOP2 and SCOP2B assignments:
            https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_scop2b_sf_uniprot.tsv.gz
            https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_scop2_uniprot.tsv.gz

        """
        urlTargetScop2 = "https://scop.mrc-lmb.cam.ac.uk/files"
        encoding = "utf-8-sig" if sys.version_info[0] > 2 else "ascii"
        fn = "scop-des-latest.txt"
        url = os.path.join(urlTargetScop2, fn)
        desL = self.__mU.doImport(url, fmt="list", uncomment=True, encoding=encoding)
        logger.info("Fetched URL is %s len %d", url, len(desL))
        #
        fn = "scop-cla-latest.txt"
        url = os.path.join(urlTargetScop2, fn)
        claL = self.__mU.doImport(url, fmt="list", uncomment=True, encoding=encoding)
        logger.info("Fetched URL is %s len %d", url, len(claL))
        #
        headerLines = self.__mU.doImport(url, fmt="list", uncomment=False, encoding=encoding)
        self.__version = headerLines[0].split(" ")[3] if headerLines else "2021-05-27"
        # JDW note cert issues with this site
        urlTargetSifts = "http://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv"
        fn = "pdb_chain_scop2b_sf_uniprot.tsv.gz"
        url = os.path.join(urlTargetSifts, fn)
        scop2bL = self.__mU.doImport(url, fmt="tdd", rowFormat="dict", uncomment=True, encoding=encoding)
        logger.info("Fetched URL is %s len %d", url, len(scop2bL))
        #
        fn = "pdb_chain_scop2_uniprot.tsv.gz"
        url = os.path.join(urlTargetSifts, fn)
        scop2L = self.__mU.doImport(url, fmt="tdd", rowFormat="dict", uncomment=True, encoding=encoding)
        logger.info("Fetched URL is %s len %d", url, len(scop2bL))
        #
        return desL, claL, scop2bL, scop2L

    def __extractNames(self, nmL):
        """ """
        rD = {}
        logger.info("Length of input name list %d", len(nmL))
        for nm in nmL:
            ff = nm.split(" ")
            rD[ff[0]] = " ".join(ff[1:])
        # self.__mU.doExport(os.path.join(self.__dirPath, "scop2-names.json"), rD, fmt="json", indent=3)
        return rD

    def __extractDomainHierarchy(self, dmL):
        """Extract the domain node identifier hierarchy from the SCOP2 representative assignment file ...

        Returns:
            dict, dict, dict, dict, dict: parent and name type dictionaries, family and superfamily assignments, and
                                          domain to superfamily mapping

            ntD[domainId] = name type TP=protein type, CL=protein class, CF=fold, SF=superfamily, FA=family
            pD[child domain identifier] = parent domain identifier
            fD[(pdbId, authAsymId)] = [(faDomId, faId, authAsymId, resBeg, resEnd),]
            sfD[(pdbId, authAsymId)] = [(sfDomId, sfId, authAsymId, resBeg, resEnd),]
            domToSfD[domSfid] = sfId

        Example assignment file:

        # SCOP release 2021-05-27
        # http://scop.mrc-lmb.cam.ac.uk
        # based on PDB release 2021-05-14
        # based on UniProt realese 2021-04-08
        # based on SIFTS release 2021-05-19
        # FA-DOMID FA-PDBID FA-PDBREG FA-UNIID FA-UNIREG SF-DOMID SF-PDBID SF-PDBREG SF-UNIID SF-UNIREG SCOPCLA
        8045703 3H8D C:1143-1264 Q64331 1143-1264 8091604 3H8D C:1143-1264 Q64331 1143-1264 TP=1,CL=1000003,CF=2001470,SF=3002524,FA=4004627
        8094330 6J56 A:1158-1282 Q9UM54 1167-1291 8094331 6J56 A:1158-1282 Q9UM54 1167-1291 TP=1,CL=1000003,CF=2001470,SF=3002524,FA=4004627
        #

        """
        # Build the parent dictionary and name node type
        ntD = {}
        pAD = {}
        pBD = {}
        pBRootD = {}
        fD = {}
        sfD = {}
        domToSfD = {}
        #
        logger.info("Length of input domain assignment list %d", len(dmL))
        for dm in dmL:
            try:
                ff = dm.split(" ")
                domFamilyId = ff[0]
                domSuperFamilyId = ff[5]
                rngL = ff[10].split(",")
                tD = {}
                for rng in rngL:
                    tL = rng.split("=")
                    tD[tL[0]] = tL[1]
                #
                # -
                # pD[tD["TP"]] = 0
                # pD[tD["CL"]] = tD["TP"]
                # pD[tD["CF"]] = tD["CL"]
                # pD[tD["SF"]] = tD["CF"]
                # pD[tD["FA"]] = tD["SF"]
                # pD[domFamilyId] = tD["FA"]
                # pD[domSuperFamilyId] = tD["SF"]
                #
                #  Represent as two trees separately rooted in protein type  and structural class
                pAD[tD["TP"]] = 0
                pAD[tD["CF"]] = tD["TP"]
                pAD[tD["SF"]] = tD["CF"]
                pAD[tD["FA"]] = tD["SF"]
                pAD[domFamilyId] = tD["FA"]
                pAD[domSuperFamilyId] = tD["SF"]
                #
                # Use this complete pBD here only for generating ID lineages, but NOT for merging with pAD
                pBD[tD["CL"]] = 0
                pBD[tD["CF"]] = tD["CL"]
                pBD[tD["SF"]] = tD["CF"]
                pBD[tD["FA"]] = tD["SF"]
                pBD[domFamilyId] = tD["FA"]
                pBD[domSuperFamilyId] = tD["SF"]
                #
                # Use pBRootD for creating tree node lists; Don't capture any lower branches to avoid re-creating redundant key:values already in pAD
                pBRootD[tD["CL"]] = 0
                pBRootD[tD["CF"]] = tD["CL"]
                #
                ntD[tD["FA"]] = "FA"
                ntD[tD["SF"]] = "SF"
                ntD[tD["CF"]] = "CF"
                ntD[tD["CL"]] = "CL"
                ntD[tD["TP"]] = "TP"
                #
                pdbId = ff[1]
                authAsymId, authSeqBeg, authSeqEnd = self.__parseAssignment(ff[2])
                if authAsymId is not None:
                    fKey = (pdbId, authAsymId)
                    fTup = (domFamilyId, tD["FA"], authAsymId, authSeqBeg, authSeqEnd)
                    if fKey not in fD:
                        fD[fKey] = []
                    if fTup not in fD[fKey]:
                        fD[fKey].append(fTup)
                pdbId = ff[6]
                authAsymId, authSeqBeg, authSeqEnd = self.__parseAssignment(ff[7])
                if authAsymId is not None:
                    sfKey = (pdbId, authAsymId)
                    sfTup = (domSuperFamilyId, tD["SF"], authAsymId, authSeqBeg, authSeqEnd)
                    if sfKey not in sfD:
                        sfD[sfKey] = []
                    if sfTup not in sfD[sfKey]:
                        sfD[sfKey].append(sfTup)
                #
                domToSfD[domSuperFamilyId] = tD["SF"]
            except Exception as e:
                logger.exception("Failing for case %r: %s", dm, str(e))
        #
        logger.info("pAD (%d) pBD (%d) pBRootD (%d) ntD (%d)", len(pAD), len(pBD), len(pBRootD), len(ntD))
        logger.info("fD (%d) sfD (%d)", len(fD), len(sfD))
        return pAD, pBD, pBRootD, ntD, fD, sfD, domToSfD

    def __parseAssignment(self, tS):
        authAsymId = authSeqBeg = authSeqEnd = None
        try:
            fL = tS.split(":")
            authAsymId = fL[0]
            rS = fL[1]
            if rS[0] == "-":
                authSeqBeg = -int(rS[1:].split("-")[0])
                authSeqEnd = int(rS[1:].split("-")[1])
            else:
                authSeqBeg = int(rS.split("-")[0])
                authSeqEnd = int(rS.split("-")[1])
        except Exception:
            pass
        return authAsymId, authSeqBeg, authSeqEnd

    def __extractScop2bSuperFamilyAssignments(self, scop2bL, domToSfD):
        """
        Extract the SCOP2B  SIFTS superfamily domain assignments for PDB structure entries.

        Returns:

         aD[(pdbId, authAsymId)] = [(sfDomId, sfId, authAsymId, resBeg, resEnd),]

        Example:

        # 2021/06/12 - 05:52 | PDB: 23.21 | UniProt: 2021.03
          PDB     CHAIN   SF_DOMID        SP_PRIMARY      RES_BEG RES_END PDB_BEG PDB_END SP_BEG  SP_END
          5id7    B       8033045 P02768  197     388     197     388     221     412
          1o9x    A       8033045 P02768  197     388     197     388     221     412
        """
        sfD = {}
        if len(scop2bL) == 0 or len(domToSfD) == 0:
            logger.error("Empty list or dict (scop2bL len %r, domToSfD len %r)", len(scop2bL), len(domToSfD))
            return sfD
        try:
            for rowD in scop2bL:
                if rowD["SF_DOMID"] in domToSfD:
                    sfD.setdefault((rowD["PDB"].upper(), rowD["CHAIN"]), []).append((rowD["SF_DOMID"], domToSfD[rowD["SF_DOMID"]], rowD["CHAIN"], rowD["PDB_BEG"], rowD["PDB_END"]))
                else:
                    logger.warning("Missing SCOP2B SF ID mapping for %r", rowD["SF_DOMID"])
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return sfD

    def __exportTreeNodeList(self, nD, pAD, pBRootD):
        """Create node list from the SCOP2 parent and name/description dictionaries.

        Exclude the root node from the tree.

        """
        #
        rootId = 0
        pL = [rootId]
        #
        logger.info("nD %d pAD %d pBRootD %d pL %r", len(nD), len(pAD), len(pBRootD), pL)
        # create child dictionary
        cD = {}
        for ctId, ptId in pAD.items():
            cD.setdefault(ptId, []).append(ctId)
        for ctId, ptId in pBRootD.items():
            cD.setdefault(ptId, []).append(ctId)
        #
        logger.debug("cD %d", len(cD))
        #
        idL = []
        for rootId in sorted(pL):
            visited = set([rootId])
            queue = collections.deque(visited)
            while queue:
                tId = queue.popleft()
                idL.append(tId)
                if tId not in cD:
                    # logger.warning("No children for scop tId %r", tId)
                    continue
                for childId in cD[tId]:
                    if childId not in visited:
                        queue.append(childId)
                        visited.add(childId)
        #
        dL = []
        for tId in idL:
            displayName = nD[tId] if tId in nD else None
            ptIdL = []
            if tId in pAD:
                ptIdL.append(pAD[tId])
            if tId in pBRootD:
                ptIdL.append(pBRootD[tId])
            lL = self.getIdLineage(tId)[1:]
            #
            # d = {'id': str(tId), 'name': displayName, 'lineage': [str(t) for t in lL], 'parents': [str(ptId)], 'depth': len(lL)}
            if tId == rootId:
                continue
            elif any([ptId == rootId for ptId in ptIdL]):
                dD = {"id": str(tId), "name": displayName, "depth": 0}
            else:
                displayName = displayName if displayName else "Domain %s" % str(tId)
                dD = {"id": str(tId), "name": displayName, "parents": ptIdL, "depth": len(lL)}
            dL.append(dD)

        return dL
