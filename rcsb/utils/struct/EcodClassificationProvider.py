##
#  File:  EcodClassificationProvider.py
#  Date:  3-Apr-2019 jdw
#
#  Updates:
#  16-Nov-2021 dwp Append additional ecod annotations for given entryId and chainId instead of overwriting
#  18-Apr-2023 aae Get version from data list directly rather than opening file twice
#  17-Mar-2026 dwp Re-work provider to handle new data format and improve robustness
#
##
"""
  Extract ECOD domain assignments, term descriptions and ECOD classification hierarchy
  from ECOD flat files.

"""

import sys
import os.path
import datetime
import logging

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
        # urlTarget = kwargs.get("ecodTargetUrl", "http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt")
        urlTarget = kwargs.get("ecodTargetUrl", "/Users/dennispiehl/Downloads/ecod.latest.domains.txt")
        # urlTarget = kwargs.get("ecodTargetUrl", "/Users/dennispiehl/rcsb/py-rcsb_utils_struct/rcsb/utils/tests-struct/sample_ecod_data.tsv")
        urlBackup = kwargs.get("ecodUrlBackupPath", "https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets/master/fall_back/ECOD/ecod.latest.domains.txt.gz")
        # self.__urlHierarchy = "http://prodata.swmed.edu/ecod/distributions/ecod.latest.hierarchy.txt"
        self.__urlHierarchy = "/Users/dennispiehl/Downloads/ecod.latest.hierarchy.txt"
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
            def dfs(node, seen):
                if node in seen:
                    return []
                seen.add(node)
                result = [node]
                for parent in self.__pD.get(node, []):
                    result.extend(dfs(parent, seen))
                return result
            #
            pList = dfs(domId, set())
            pList.reverse()
        except Exception as e:
            logger.exception("Failing for %r with %s", domId, str(e))
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
        return self.__exportTreeNodeList(self.__pD, self.__nD)

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
            tS = datetime.datetime.now().isoformat()  # Or, grab from domains flat file header?
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
            ok = fU.get(urlTarget, fp)
            if not ok:
                return None
        #
        nmdLUnfiltered = self.__mU.doImport(fp, fmt="tdd", rowFormat='dict')
        keys_to_keep = {
            "ecod_domain_id",
            "f_id",
            "pdb",
            "chain",
            "pdb_range",
            "architecture_name",
            "x_name",
            "h_name",
            "t_name",
            "f_name",
        }

        nmL = [{k: d[k].replace('"', "").strip() for k in keys_to_keep if k in d} for d in nmdLUnfiltered]

        with open(fp, encoding="utf-8") as f:
            version = next(
                (line.partition(":")[2].strip() for line in f if line.startswith("# Version:")),
                None,
            )
        if version is None:
            logger.error("No version line found in file %r", fp)
        self.__version = version
        logger.info("Parsed version from header %r", version)

        fU.remove(fp)
        return nmL

    def __extractDomainHierarchy(self, nmL):
        """
        # ECOD Domain List
        # Version: develop294
        # Generated: 2026-02-22 20:10:01
        #
        uid	ecod_domain_id	manual_rep	f_id	pdb	chain	pdb_range	seqid_range	architecture_name	x_name	h_name	t_name	f_name	assembly_id	domain_id_short	range_count	arch_manual	x_manual	h_manual	t_manual	f_manual	valid_structure	ligand_binding
        0	e2nmzA1	TRUE	1.1.1.3	2nmz	A	A:1-99	A:1-99	beta barrels	cradle loop barrel	RIFT-related	acid protease	RVP				FALSE	FALSE	FALSE	FALSE	TRUE	TRUE	FALSE
        1	e1hvcA1	FALSE	1.1.1.3	1hvc	A	A:1B-99A	A:1-203	beta barrels	cradle loop barrel	RIFT-related	acid protease	RVP				FALSE	FALSE	FALSE	FALSE	FALSE	TRUE	FALSE
        2	e4fivA1	TRUE	1.1.1.3	4fiv	A	A:4-116	A:1-113	beta barrels	cradle loop barrel	RIFT-related	acid protease	RVP				FALSE	FALSE	FALSE	FALSE	TRUE	TRUE	FALSE
        5084102	UPI003767874A_nD1	FALSE	2003.1.1.20				1-335	a/b three-layered sandwiches	Rossmann-like	Rossmann-related	NAD(P)-binding Rossmann-fold domains	Epimerase				FALSE	FALSE	FALSE	FALSE	FALSE	TRUE	FALSE
        5084103	UPI00376788AD_nD1	FALSE	11.1.1.0				1-100	beta sandwiches	Immunoglobulin-like beta-sandwich	Immunoglobulin-related	Immunoglobulin/Fibronectin type III/E set domains/PapD-like					FALSE	FALSE	FALSE	FALSE	FALSE	TRUE	FALSE
        4891150	e9yiuG1	False	315.1.1.1	9yiu	G	G:1-114	G:1-114		Tautomerase/MIF-like	Tautomerase/MIF	Tautomerase/MIF	MIF				False	False	False	False	False	True	False
        1234748	e3id6A3	False	606.1.1.1	3id6	A	A:133-259	A:133-259	alpha complex topology	Nop N-terminal domain	Nop N-terminal domain	Nop N-terminal domain	Nop				False	False	False	False	False	True	False

        Args:
            nmL (list): list of dictionaries containing all domain data (filtered to only include the used portions)

        Returns:
            pD (dict):   dictionary containing all children IDs as keys and all possible parent IDs as values
            nmD (dict):  dictionary mapping of ID (key) to full name (including leading hierarchy letter)
            ntD (dict):  dictionary mapping of ID (key) to hierarchy letter (value)
            pdbD (dict): dictionary containing all PDB ID & chain combos as keys and the list of associated ECOD assignment & alignment info tuples as values

        Example output:
            # pD:   {"102.1.1.146": ['3407.1.1'], "304.4.1.66": ['304.7.1'], ...}
            # nmD:  {'a': 'A: special', 'a.1': 'A: beta barrels', 'a.10': 'A: alpha complex topology', 'a.11': 'A: a+b two layers', ...}
            # ntD:  {"a": "A", "a.1": "A", "a.10": "A", "a.11": "A"}
            # pdbD:
                    ('2nmz', 'A') [('e2nmzA1', '1.1.1.3', 'A', 1, 99)]
                    ('1hvc', 'A') []                                ## IS THIS OK?????
                    ('4fiv', 'A') [('e4fivA1', '1.1.1.3', 'A', 4, 116)]
        """
        pD = {}
        parentChildTupleSet = set()
        pdbD = {}
        #
        # First read hierarchy definition file
        ntD, nmD, letterNameToIdD = self.__parseHierarchy()
        logger.info("Length of ntD %r, nmD %r, letterNameToIdD %r", len(ntD), len(nmD), len(letterNameToIdD))
        #TODO: Why are nmD and letterNameToIdD different? they should be mirrors? -- Length of ntD 47986, nmD 47986, letterNameToIdD 37984
        #
        # Now parse the data file using the hierarchy id:name mappings
        logger.info("Length of input ECOD name list %d", len(nmL))
        for nm in nmL:
            # ecodId is the linkable identifier -
            ecodId = nm["ecod_domain_id"]
            f_id = nm["f_id"]
            entryId = nm.get("pdb").lower()
            if not entryId:
                continue
            authAsymId = nm["chain"]
            resRange = nm["pdb_range"]  # TO REVIEW: can we use seqid_range instead or in addition?
            #
            # TODO: Compare if below comment still holds or not, and whether I need to up-populate other fields
            #  There are no unique identifiers published for the internal elements of the hierarchy
            #   so these are assigned here similar to scop - There are also many unnamed nodes
            #   that are conventionally filled in from the leaf levels of the tree...
            #  {"A": "Architecture", "X": "Possible Homology", "H": "Homology", "T": "Topology", "F": "Family"}

            aN = nm.get("architecture_name")
            xN = nm.get("x_name")
            hN = nm.get("h_name")
            tN = nm.get("t_name")
            fN = nm.get("f_name")
            if not aN:
                logger.debug("ecodId %r entryId %r has no architecture group - skipping", ecodId, entryId)
                continue

            if fN:


            aId = letterNameToIdD.get(f"A: {aN}")
            xId = letterNameToIdD.get(f"X: {xN}") if xN else None
            hId = letterNameToIdD.get(f"H: {hN}") if hN else None
            tId = letterNameToIdD.get(f"T: {tN}") if tN else None
            fId = letterNameToIdD.get(f"F: {fN}") if fN else None

            if aId and xId:
                parentChildTupleSet.add((aId, xId))
            if xId and hId:
                parentChildTupleSet.add((xId, hId))
            if hId and tId:
                parentChildTupleSet.add((hId, tId))
            if tId and fId:
                parentChildTupleSet.add((tId, fId))
            #
            rL = self.__parseRanges(resRange)
            if (entryId, authAsymId) not in pdbD:
                pdbD[(entryId, authAsymId)] = [(ecodId, fId, t[0], t[1], t[2]) for t in rL]
            else:
                for t in rL:
                    pdbD[(entryId, authAsymId)].append((ecodId, fId, t[0], t[1], t[2]))

        parentChildTupleList = list(parentChildTupleSet)

        # create a dictionary to store the parents of each child
        for parent, child in parentChildTupleList:
            if child not in pD:
                pD[child] = []
            if parent not in pD[child]:
                pD[child].append(parent)

        return pD, nmD, ntD, pdbD

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

    def __exportTreeNodeList(self, childToParentD, idNameMapD):
        """Create tree node list in the format of:

        {'id': '3193', 'name': 'X: HopAB effectors Pto-binding domain-related', 'parents': ['a.6']}
        {'id': '2002.1.1.190', 'name': 'F: DUF5696', 'parents': ['2002.1.1']}
        {'id': '2484.3.1.7', 'name': 'F: Peptidase_M24', 'parents': ['283.1.1']}
        {'id': '2484.1.1.96', 'name': 'F: DNA_polI_exo1', 'parents': ['2484.1.1']}
        """
        #
        print("IN EXPORT TREE")
        dL = []
        for child, parentL in childToParentD.items():
            if parentL:
                tD = {"id": child, "name": idNameMapD[child], "parents": parentL}
            else:
                tD = {"id": child, "name": idNameMapD[child]}
            dL.append(tD)
        # x=0
        # for d in dL:
        #     print(d)
        #     x+=1
        #     if x == 100:
        #         break

        return dL

    # def __exportTreeNodeList(self, pD):
    #     """Create node list from name dictionary and lineage dictionaries."""
    #     #
    #     rootId = 0
    #     pL = [rootId]
    #     #
    #     logger.info("pD %d pL %r", len(pD), pL)
    #     # --
    #     #
    #     # create child dictionary
    #     cD = {}
    #     for ctId, ptId in pD.items():
    #         cD.setdefault(ptId, []).append(ctId)
    #     #
    #     logger.info("cD %d", len(cD))
    #     #
    #     idL = []
    #     for rootId in sorted(pL):
    #         visited = set([rootId])
    #         queue = collections.deque(visited)
    #         while queue:
    #             tId = queue.popleft()
    #             idL.append(tId)
    #             if tId not in cD:
    #                 # logger.debug("No children for Ecod tId %s", tId)
    #                 continue
    #             for childId in cD[tId]:
    #                 if childId not in visited:
    #                     queue.append(childId)
    #                     visited.add(childId)
    #     #
    #     dL = []
    #     for tId in idL:
    #         displayName = self.getName(tId)
    #         ptId = pD[tId] if tId in pD else None
    #         lL = self.getIdLineage(tId)[1:]
    #         #
    #         if tId == rootId:
    #             continue
    #         elif ptId == rootId:
    #             dD = {"id": str(tId), "name": displayName, "depth": 0}
    #         else:
    #             dD = {"id": str(tId), "name": displayName, "parents": [str(ptId)], "depth": len(lL)}
    #         dL.append(dD)

    #     return dL

    def __parseHierarchy(self):
        """
        Parse the ECOD hierarchy definition file and generate mapping dictionaries
        for level IDs, names, and hierarchy letters.

        Reads the hierarchy file from `self.__urlHierarchy` using `self.__mU.doImport`
        and returns four dictionaries for convenient lookups.

        Returns:
            idToLetterD (dict[str, str]):
                Maps each level ID to its hierarchy letter (A, X, H, T, F).
                Example: {'a.1': 'A', 'a.10': 'A'}

            idToLetterNameD (dict[str, str]):
                Maps each level ID to the combined letter + name string.
                Example: {'a.1': 'A: beta barrels'}

            letterNameToIdD (dict[str, str]):
                Maps the combined letter + name string back to its level ID.
                Example: {'A: beta barrels': 'a.1'}

        Raises:
            ValueError: If any line in the hierarchy file has fewer than 3 columns.

        Notes:
            - Lines with missing names (except Family "F" levels) are labeled as "special".
            - The hierarchy file is expected to be tab-delimited, with at least:
                [level_letter, level_id, level_name]
        """
        idToLetterD = {}  # ntD
        idToLetterNameD = {}  # nmD
        letterNameToIdD = {}
        hlL = self.__mU.doImport(self.__urlHierarchy, fmt="tdd", rowFormat='list')
        # ['A', 'a', '', '', '7241']
        # ['A', 'a.1', 'beta barrels', '', '209685']

        for hli in hlL:

            # Defensive: ensure at least 3 columns
            if len(hli) < 3:
                raise ValueError(f"Malformed line: {hli}")

            level_letter = hli[0].strip()
            level_id = hli[1].strip()
            level_name = hli[2].strip()
            if not level_name:  # should we do this?
                continue
            # if level_letter != "F" and level_name == "":
            #     print(level_letter, level_id)
            #     level_name = "special"

            if level_id not in idToLetterNameD:
                idToLetterD[level_id] = level_letter
                idToLetterNameD[level_id] = f"{level_letter}: {level_name}"
            if f"{level_letter}: {level_name}" not in letterNameToIdD:
                letterNameToIdD[f"{level_letter}: {level_name}"] = level_id

        # print("analyzing mirrored data")
        # for k,v in idToLetterNameD.items():
        #     for k2,v2 in idToLetterNameD.items():
        #         if k != k2 and v == v2:
        #             print(k, k2, v, v2)
        # LOTS OF CASES:
            # 3081 3401 X: Envelope glycoprotein GP1 X: Envelope glycoprotein GP1
            # 3401 3081 X: Envelope glycoprotein GP1 X: Envelope glycoprotein GP1
            # 3407 606 X: Nop N-terminal domain X: Nop N-terminal domain
            # 606 3407 X: Nop N-terminal domain X: Nop N-terminal domain
            # 327.1.1.3 101.43.1.2 F: Transposase_mut F: Transposase_mut
            # 327.1.1.3 2484.1.1.199 F: Transposase_mut F: Transposase_mut
            # 327.1.1.3 4120.1.1.41 F: Transposase_mut F: Transposase_mut
        return idToLetterD, idToLetterNameD, letterNameToIdD


# Every line in the file has a 'f_id', which corresponds to the smallest labeled hierarchy. Will need to use that to work backwards to form lineage, instead of the Name:ID mapping
# Seems to be possible to just chop away at the last decimal, at least up until A, which uses its own ID....

# (.venv) dennis@DennisMBP [26-03-17 18:11:52] ~/rcsb/py-rcsb_utils_struct/rcsb/utils/tests-struct % grep 'Nop N-terminal domain' '/Users/dennispiehl/Downloads/ecod.latest.domains.txt' | cut -f4,9,10,11,12,13 | sort -u
# 3407.1.1.0        mixed a+b and a/b       Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain
# 3407.1.1.1        Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   NOP5NT
# 3407.1.1.1        mixed a+b and a/b       Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   NOP5NT
# 3407.1.1.2        mixed a+b and a/b       Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   Nop5_56-rel_N_Arc
# 3407.1.1.3        mixed a+b and a/b       Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   Nop5_N
# 3407.1.1.4        Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   NOP5NT, Nop5_56-rel_N_Arc
# 3407.1.1.4        mixed a+b and a/b       Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   NOP5NT, Nop5_56-rel_N_Arc
# 3407.1.1.5        mixed a+b and a/b       Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   Nop
# 606.1.1.0         Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain
# 606.1.1.0         alpha complex topology  Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain
# 606.1.1.1         Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   Nop
# 606.1.1.1         alpha complex topology  Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   Nop
# 606.1.1.10        Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   DOG1
# 606.1.1.12        Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   Rx_N
# 606.1.1.13        Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   DUF3475
# 606.1.1.14        Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   Acyl_transf_3
# 606.1.1.15        Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   NFACT-C
# 606.1.1.16        Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   Transposase_20
# 606.1.1.16        alpha complex topology  Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   Transposase_20
# 606.1.1.18        alpha complex topology  Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   Cys_rich_FGFR
# 606.1.1.3         alpha complex topology  Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   DEDD_Tnp_IS110
# 606.1.1.7         Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   TYA
# 606.1.1.8         Nop N-terminal domain   Nop N-terminal domain   Nop N-terminal domain   APO_RNA-bind

# pD:   {100001: 0, 200001: 100001, 300001: 200001, 400001: 300001, 500001: 400001, 100002: 0, 200002: 100002, 300002: 200002, 400002: 300002, 500002: 400002, 100003: 0, 200003: 100003, 300003: 200003, 400003: 300005, 500003: 400003, 100004: 0, 200004: 100004, 300004: 200004, 100005: 0, 200005: 100005, 300005: 200005}

# nmD:   {100001: 'A: beta barrels', 200001: 'X: cradle loop barrel', 300001: 'H: RIFT-related', 400001: 'T: acid protease', 500001: 'F: RVP', 100002: 'A: a/b barrels', 200002: 'X: Ten stranded beta/alpha barrel', 300002: 'H: Ten stranded beta/alpha barrel', 400002: 'T: Ten stranded beta/alpha barrel', 500002: 'F: DUF711', 100003: 'A: alpha complex topology', 200003: 'X: Nop N-terminal domain', 300003: 'H: Nop N-terminal domain', 400003: 'T: Nop N-terminal domain', 500003: 'F: Nop', 100004: 'A: alpha complex topology0000', 200004: 'X: Nop N-terminal domain', 300004: 'H: Nop N-terminal domain', 100005: 'A: alpha complex topology000', 200005: 'X: Nop N-terminal domain', 300005: 'H: Nop N-terminal domain'}

# ntD:   {100001: 'A', 200001: 'X', 300001: 'H', 400001: 'T', 500001: 'F', 100002: 'A', 200002: 'X', 300002: 'H', 400002: 'T', 500002: 'F', 100003: 'A', 200003: 'X', 300003: 'H', 400003: 'T', 500003: 'F', 100004: 'A', 200004: 'X', 300004: 'H', 100005: 'A', 200005: 'X', 300005: 'H'}

# pdbD:   {('2nmz', 'A'): [('e2nmzA1', 500001, 'A', 1, 99)], ('1hvc', 'A'): [], ('4fiv', 'A'): [('e4fivA1', 500001, 'A', 4, 116)], ('2rsp', 'A'): [('e2rspA1', 500001, 'A', 1, 124)], ('2fmb', 'A'): [('e2fmbA1', 500001, 'A', 1, 104)], ('2ha9', 'A'): [('e2ha9A1', 500002, 'A', 1, 440)], ('2ha9', 'B'): [('e2ha9B1', 500002, 'B', 1, 440)], ('3id6', 'A'): [('e3id6A3', 500003, 'A', 133, 259), ('e3id6A2', 500003, 'A', 133, 259), ('e3id6A1', 500003, 'A', 133, 259)]}


# pD:   {100001: 0, 200001: 100001, 300001: 200001, 400001: 300001, 500001: 400001, 100002: 0, 200002: 100002, 300002: 200002, 400002: 300002, 500002: 400002, 100003: 0, 200003: 100003, 300003: 200003, 400003: 300005, 500003: 400003, 100004: 0, 200004: 100004, 300004: 200004, 100005: 0, 200005: 100005, 300005: 200005}
# nmD:   {100001: 'A: beta barrels', 200001: 'X: cradle loop barrel', 300001: 'H: RIFT-related', 400001: 'T: acid protease', 500001: 'F: RVP', 100002: 'A: a/b barrels', 200002: 'X: Ten stranded beta/alpha barrel', 300002: 'H: Ten stranded beta/alpha barrel', 400002: 'T: Ten stranded beta/alpha barrel', 500002: 'F: DUF711', 100003: 'A: alpha complex topology', 200003: 'X: Nop N-terminal domain', 300003: 'H: Nop N-terminal domain', 400003: 'T: Nop N-terminal domain', 500003: 'F: Nop', 100004: 'A: alpha complex topology0000', 200004: 'X: Nop N-terminal domain', 300004: 'H: Nop N-terminal domain', 100005: 'A: alpha complex topology000', 200005: 'X: Nop N-terminal domain', 300005: 'H: Nop N-terminal domain'}
# ntD:   {100001: 'A', 200001: 'X', 300001: 'H', 400001: 'T', 500001: 'F', 100002: 'A', 200002: 'X', 300002: 'H', 400002: 'T', 500002: 'F', 100003: 'A', 200003: 'X', 300003: 'H', 400003: 'T', 500003: 'F', 100004: 'A', 200004: 'X', 300004: 'H', 100005: 'A', 200005: 'X', 300005: 'H'}
# pdbD:   {('2nmz', 'A'): [('e2nmzA1', 500001, 'A', 1, 99)], ('1hvc', 'A'): [], ('4fiv', 'A'): [('e4fivA1', 500001, 'A', 4, 116)], ('2rsp', 'A'): [('e2rspA1', 500001, 'A', 1, 124)], ('2fmb', 'A'): [('e2fmbA1', 500001, 'A', 1, 104)], ('2ha9', 'A'): [('e2ha9A1', 500002, 'A', 1, 440)], ('2ha9', 'B'): [('e2ha9B1', 500002, 'B', 1, 440)], ('3id6', 'A'): [('e3id6A3', 500003, 'A', 133, 259), ('e3id6A2', 500003, 'A', 133, 259), ('e3id6A1', 500003, 'A', 133, 259)]}



## FRom production before new ECOD data:
## Note: Use of "F: ...", .and "annotation_id" is the ecod_domain_id
# {
#   "annotation_id": "e2ha9B1",
#   "assignment_version": "1.6",
#   "description": null,
#   "name": "DUF711",
#   "provenance_source": "ECOD",
#   "annotation_lineage": [
#     {
#       "id": "100016",
#       "name": "A: a/b barrels",
#       "depth": 1
#     },
#     {
#       "id": "201880",
#       "name": "X: Ten stranded beta/alpha barrel (From Topology)",
#       "depth": 2
#     },
#     {
#       "id": "303015",
#       "name": "H: Ten stranded beta/alpha barrel (From Topology)",
#       "depth": 3
#     },
#     {
#       "id": "403171",
#       "name": "T: Ten stranded beta/alpha barrel",
#       "depth": 4
#     },
#     {
#       "id": "511542",
#       "name": "F: DUF711",
#       "depth": 5
#     }
#   ]
