##
#  File:           CathDictionary.py
#  Original Date:  12-Nov-2007 jdw
# 
#  Updated:
#
#  26-Apr-2009 jdw  Add additional methods to extract classification hierarchy.
#   9-May-2011 jdw  import to sbkb framework.
##
"""
  Extract CATH term descriptions and CATH classification hierarchy
  from CATH flat files.
  
"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.001"

import sys

from sbkb.utils.ConfigInfo                 import ConfigInfo
from sbkb.fetchers.resource.ResourceConfig import ResourceConfig

#
class CathDictionary:
    """Generate a dictionary of CATH classification terms.

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
        self.descrPath=self.__getDescrFilePath()
        self.__cathTermDict={}
        #
        self.__makeCathTermDict()
        #

    def __getDescrFilePath(self):
        """ Get the path to CATH description file CathDomainDescriptionFile.v3.2.0
        """
        fid=self.__cI.get('CATH_CLASS_DESCRIPTIONS_FID')
        fPath=self.__rC.getLocalPath(fid[0],fid[1])
        if (not self.__rC.localPathExists(fid[0],fid[1])):
            self.__lfh.write("+CathDictionary(__getDescrFilePath) CATH class description file not found %s\n" % fPath)
        return (fPath)

    def getCathTermDict(self):
        return (self.__cathTermDict)

    def getCathText(self,domain):
        if (self.__cathTermDict.has_key(domain.upper())):
            d=self.__cathTermDict[domain.upper()]
            textC  = "Class: " + d['CLASS']
            textC += " |Architecture: " + d['ARCH']
            textC += " |Topology: " + d['TOPOL']
            textC += " |Homology: " + d['HOMOL']                         
            return textC
        else:
            return None

    def getCathTextTdd(self,domain):
        if (self.__cathTermDict.has_key(domain.upper())):
            d=self.__cathTermDict[domain.upper()]
            textC  = d['CLASS']
            textC += '\t'
            textC += d['ARCH']
            textC += '\t'            
            textC += d['TOPOL']
            textC += '\t'            
            textC += d['HOMOL']                         
            return textC
        else:
            return None

    def getCathCode(self,domain):
        if (self.__cathTermDict.has_key(domain.upper())):
            d=self.__cathTermDict[domain.upper()]
            return d['CATHCODE']
        else:
            return None

    def getCathHomology(self,domain):
        if (self.__cathTermDict.has_key(domain.upper())):
            d=self.__cathTermDict[domain.upper()]
            textC  = "Class: " + d['CLASS']+ " |Architecture: " + d['ARCH'] + " |Topology: " + d['TOPOL'] + " |Homology: " + d['HOMOL']
            tCode = str(d['CATHCODE'])
            return (tCode,textC)
        else:
            return (None,None)

    def getCathTopology(self,domain):
        if (self.__cathTermDict.has_key(domain.upper())):
            d=self.__cathTermDict[domain.upper()]
            textC  = "Class: " + d['CLASS']+ " |Architecture: " + d['ARCH'] + " |Topology: " + d['TOPOL']
            cCode=str(d['CATHCODE']).split('.')
            tCode=".".join(cCode[:3])
            return (tCode,textC)
        else:
            return (None,None)

    def getCathArchitecture(self,domain):
        if (self.__cathTermDict.has_key(domain.upper())):
            d=self.__cathTermDict[domain.upper()]
            textC  = "Class: " + d['CLASS']+ " |Architecture: " + d['ARCH']
            cCode=str(d['CATHCODE']).split('.')
            tCode=".".join(cCode[:2])
            return (tCode,textC)
        else:
            return (None,None)

    def getCathClass(self,domain):
        if (self.__cathTermDict.has_key(domain.upper())):
            d=self.__cathTermDict[domain.upper()]
            textC  = "Class: " + d['CLASS']
            cCode=str(d['CATHCODE']).split('.')
            tCode=str(cCode[0])
            return (tCode,textC)
        else:
            return (None,None)

        
    def __makeCathTermDict(self):
        """Return a dictionary containing CATH terms - 
        """
        self.__cathTermDict={}
        entry={}        
        f=open(self.descrPath,'r')
        isFirst = True
        for line in f.xreadlines():
            if line.startswith("#"):
                continue
            f1=str(line[:9]).strip()
            f2=str(line[10:-1]).strip()
            if (f1 == "//"):
                # new entry
                if isFirst:
                    isFirst = False
                else:
                    self.__cathTermDict[entry['DOMAIN'].upper()]= entry                
                entry={}
                entry['DOMAIN']=None                
                entry['CATHCODE']=None
                entry['ARCH']=None
                entry['CLASS']=None
                entry['TOPOL']=None
                entry['HOMOL']=None
            else:
                if (f1 == "DOMAIN"):
                    entry['DOMAIN']=f2                
                elif (f1 == "CATHCODE"):
                    entry["CATHCODE"]=f2
                elif (f1 == "CLASS"):
                    entry['CLASS']=f2
                elif (f1 == "ARCH"):
                    entry['ARCH']=f2                    
                elif (f1 == "TOPOL"):
                    entry['TOPOL']=f2                                        
                elif (f1 == "HOMOL"):                                        
                    entry['HOMOL']=f2                    
                    
            
        if (self.__verbose): self.__lfh.write("CATH term dictionary length %d \n" % len(self.__cathTermDict) )
        return self.__cathTermDict

    
if __name__ == '__main__':
    pass
#
