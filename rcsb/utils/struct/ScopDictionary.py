##
#  File:           ScopDictionary.py
#  Original Date:  12-Nov-2007 JDW
#
#
#
#  Updates:
#  26-Apr-2009 jdw  Add methods to extract class, fold, superfamily and family
#                   details.
#   9-May-2011 jdw  incorporate within sbkb framework 
##
"""
  Extract SCOP term descriptions from scop flat files and
  the SCOP classification hierarchy.
"""

__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.001"

import sys
from copy import deepcopy

from sbkb.utils.ConfigInfo                 import ConfigInfo
from sbkb.fetchers.resource.ResourceConfig import ResourceConfig

class ScopDictionary:
    """Generate a dictionary of SCOP classification terms. 
       
    """
    def __init__(self,configLevel='DEV',verbose=True,log=sys.stderr):
        #
        self.__verbose=verbose
        self.__lfh=log
        self.__configLevel=configLevel
        #
        self.__cI=ConfigInfo(self.__configLevel)
        self.__rC=ResourceConfig(configLevel=self.__configLevel,verbose=self.__verbose,log=self.__lfh)
    
        #
        self.__classFilePath=self.__getAssignFilePath()
        self.__descrFilePath=self.__getDescrFilePath()

        self.__scopTermDict={}
        self.__scopDict={}        
        #
        self.__makeScopTermDict()
        self.__makeScopDict()        
        #
        #

    def __getDescrFilePath(self):
        """ Get the path to SCOP class description file. 
        """
        fid=self.__cI.get('SCOP_CLASS_DESCRIPTIONS_FID')
        fPath=self.__rC.getLocalPath(fid[0],fid[1])
        if (not self.__rC.localPathExists(fid[0],fid[1])):
            self.__lfh.write("+ScopDictionary(__getDescrFilePath) SCOP class description file not found %s\n" % fPath)
        return (fPath)

    def __getAssignFilePath(self):
        """ Get the path to SCOP class term definitions 
        """
        fid=self.__cI.get('SCOP_CLASS_ASSIGNMENTS_FID')
        fPath=self.__rC.getLocalPath(fid[0],fid[1])
        if (not self.__rC.localPathExists(fid[0],fid[1])):
            self.__lfh.write("+ScopDictionary(__getClassFilePath) SCOP class definition file not found %s\n" % fPath)
        return (fPath)


    def getScopTermDict(self):
        return (self.__scopTermDict)

    def getScopDict(self):
        return (self.__scopDict)    

    def getTermDescription(self,scopAccession):
        if (self.__scopTermDict.has_key(scopAccession)):
            d= self.__scopTermDict[scopAccession]
            return d['DESCRIPTION']
        else:
            return None

    def getTermClass(self,scopAccession):
        if (self.__scopTermDict.has_key(scopAccession)):
            d= self.__scopTermDict[scopAccession]
            return d['CLASS']
        else:
            return None

    def getScopText(self,scopAccession):
        if (self.__scopDict.has_key(scopAccession)):
            d=self.__scopDict[scopAccession]
            return d['SCOP_TEXT']
        else:
            return None
        
    def __makeScopTermDict(self):
        """Return a dictionary containing SCOP terms - 
        """
        self.__scopTermDict={}
        f=open(self.__descrFilePath,'r')
        for line in f:
            if line.startswith("#"):
                continue
            fieldList=line.split('\t')
            entry = {}
            entry['ACCESSION']    = str(fieldList[0]).strip()
            entry['TYPE']         = str(fieldList[1]).strip()
            entry['CLASS']        = str(fieldList[2]).strip()
            entry['OLD_CODE']     = str(fieldList[3]).strip()
            entry['DESCRIPTION']  = str(fieldList[4]).strip()
            self.__scopTermDict[entry['ACCESSION']]= entry
        #sys.stderr.write("SCOP term dictionary length %d \n" % len(self.__scopTermDict) )
        f.close()
        return self.__scopTermDict

    def __makeScopDict(self):
        """Return a dictionary containing SCOP annotations - 
        """
        self.__scopDict={}
        f=open(self.__classFilePath,'r')
        for line in f:
            if line.startswith("#"):
                continue
            fieldList=line.split('\t')
            entry = {}
            entry['OLD_CODE']     = str(fieldList[0]).strip()
            entry['PDB_ID']       = str(fieldList[1]).strip()
            entry['CLASS']        = str(fieldList[3]).strip()
            entry['ACCESSION']    = str(fieldList[4]).strip()

            scopAnn = str(fieldList[5]).strip()
            #sys.stderr.write("%s\n" % scopAnn)            
            #
            td= {}
            vvList = scopAnn.split(",")
            if len(vvList) > 0:
                for vv in vvList:
                    ll = vv.split('=')
                    tk = str(ll[0]).strip()
                    tv = str(ll[1]).strip()
                    td[tk] = tv
                    #sys.stderr.write("%s %s\n" % (tk,tv))                               
            else:
                sys.stderr.write("Error for %s - %s\n" % (entry['OLD_CODE'],scopAnn))

            scopText=""
            if td.has_key("cl"):
                entry['CLASS']= self.getTermDescription(td["cl"])
                scopText+="Class: %s" % self.getTermDescription(td["cl"])
                entry['CLASS_TEXT']=scopText
                entry['CLASS_CLASS']=self.getTermClass(td["cl"])
                entry['CLASS_SID']=td["cl"]                
            else:
                entry['CLASS']=None
                entry['CLASS_TEXT']=None
                entry['CLASS_CLASS']=None
                entry['CLASS_SID']=None                                                        

            if td.has_key("cf"):
                entry['FOLD']=self.getTermDescription(td["cf"])
                scopText+=" |Fold: %s" % self.getTermDescription(td["cf"])
                entry['FOLD_TEXT']=scopText
                entry['FOLD_CLASS']=self.getTermClass(td["cf"])
                entry['FOLD_SID']=td["cf"]                                
            else:
                entry['FOLD']=None
                entry['FOLD_TEXT']=None
                entry['FOLD_CLASS']=None
                entry['FOLD_SID']=None                                                
                
            if td.has_key("sf"):
                entry['SUPER_FAMILY']= self.getTermDescription(td["sf"])                                
                scopText+=" |Super Family: %s" % self.getTermDescription(td["sf"])
                entry['SUPER_FAMILY_TEXT']=scopText
                entry['SUPER_FAMILY_CLASS']=self.getTermClass(td["sf"])
                entry['SUPER_FAMILY_SID']=td["sf"]                
            else:
                entry['SUPER_FAMILY']=None
                entry['SUPER_FAMILY_TEXT']=None
                entry['SUPER_FAMILY_CLASS']=None
                entry['SUPER_FAMILY_SID']=None                                                
                
            if td.has_key("fa"):
                entry['FAMILY']=self.getTermDescription(td["fa"])                      
                scopText+=" |Family: %s" % self.getTermDescription(td["fa"])
                entry['FAMILY_TEXT']=scopText
                entry['FAMILY_CLASS']=self.getTermClass(td["fa"])
                entry['FAMILY_SID']=td["fa"]                                                                
            else:
                entry['FAMILY']=None
                entry['FAMILY_TEXT']=None
                entry['FAMILY_CLASS']=None
                entry['FAMILY_SID']=None                                                

            if td.has_key("dm"):
                entry['DOMAIN']=self.getTermDescription(td["dm"])
                scopText+=" |Domain: %s" % self.getTermDescription(td["dm"])
                entry['DOMAIN_TEXT']=scopText
                entry['DOMAIN_CLASS']=self.getTermClass(td["dm"])
                entry['DOMAIN_SID']=td["dm"]                   
            else:
                entry['DOMAIN']=None
                entry['DOMAIN_TEXT']=None
                entry['DOMAIN_CLASS']=None
                entry['DOMAIN_SID']=None                                                

            if td.has_key("sp"):
                entry['SPECIES']=self.getTermDescription(td["sp"])                                        
                scopText+="|Species: %s" % self.getTermDescription(td["sp"])                                        
            else:
                entry['SPECIES']=None

            entry['SCOP_TEXT']=scopText

            chnRes = str(fieldList[2]).strip()
            #sys.stderr.write("chnRes %s\n" % chnRes)
            for cR in chnRes.split(','):
                tEntry=deepcopy(entry)
                tEntry['PDB_CHAIN_ID']=None
                tEntry['PDB_BEG_RES'] =None
                tEntry['PDB_END_RES'] =None
                if (cR.find(':') > 0):
                    crL=cR.split(':')
                    tEntry['PDB_CHAIN_ID']=crL[0]
                    if (len(crL) > 0 and len(crL[1]) > 0):
                        resL=crL[1].split('-')
                        if len(resL[0]) > 0:
                            tEntry['PDB_BEG_RES']=resL[0]
                        if len(resL) > 0 :
                            tEntry['PDB_END_RES']=resL[1]
                else:
                    tEntry['PDB_CHAIN_ID']=crL                    
                self.__scopDict[entry['ACCESSION']]= tEntry
        #        sys.stderr.write("SCOP dictionary length %d \n" % len(self.__scopDict) )
        return self.__scopDict


    def getClassDescription(self,scopAccession):
        if (self.__scopDict.has_key(scopAccession)):
            return (self.__scopDict[scopAccession]['CLASS_SID'],self.__scopDict[scopAccession]['CLASS_CLASS'],self.__scopDict[scopAccession]['CLASS_TEXT'])
            #
        else:
            return (None,None)

    def getFoldDescription(self,scopAccession):
        if (self.__scopDict.has_key(scopAccession)):
            return (self.__scopDict[scopAccession]['FOLD_SID'],self.__scopDict[scopAccession]['FOLD_CLASS'],self.__scopDict[scopAccession]['FOLD_TEXT'])
            #
        else:
            return (None,None)

    def getSuperFamilyDescription(self,scopAccession):
        if (self.__scopDict.has_key(scopAccession)):
            return (self.__scopDict[scopAccession]['SUPER_FAMILY_SID'],self.__scopDict[scopAccession]['SUPER_FAMILY_CLASS'],self.__scopDict[scopAccession]['SUPER_FAMILY_TEXT'])
            #
        else:
            return (None,None)

    def getFamilyDescription(self,scopAccession):
        if (self.__scopDict.has_key(scopAccession)):
            return (self.__scopDict[scopAccession]['FAMILY_SID'],self.__scopDict[scopAccession]['FAMILY_CLASS'],self.__scopDict[scopAccession]['FAMILY_TEXT'])
            #
        else:
            return (None,None)


    def getDomainDescription(self,scopAccession):
        if (self.__scopDict.has_key(scopAccession)):
            return (self.__scopDict[scopAccession]['DOMAIN_SID'],self.__scopDict[scopAccession]['DOMAIN_CLASS'],self.__scopDict[scopAccession]['DOMAIN_TEXT'])
            #
        else:
            return (None,None)



if __name__ == '__main__':
    pass
#
