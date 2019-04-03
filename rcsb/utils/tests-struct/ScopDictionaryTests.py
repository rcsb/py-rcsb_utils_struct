##
# File:    ScopDictionaryTests.py
# Author:  jdw
# Date:    9-May-2011
# Version: 0.001
#
# Updates:
#
##
"""
Test cases for operations that read SCOP term and class data from flat files -

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.001"

import sys, unittest, traceback
import sys, time, os, os.path, shutil


from sbkb.utils.ConfigInfo             import ConfigInfo
from sbkb.readers.scop.ScopDictionary  import ScopDictionary

class ScopDictionaryTests(unittest.TestCase):
    def setUp(self):
        self.__verbose=True
        self.__lfh=sys.stderr

    def tearDown(self):
        pass

    def testGetScopTermDictionary(self): 
        """ Get dictionary of SCOP terms 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            cD=ScopDictionary(configLevel='DEV',verbose=self.__verbose,log=self.__lfh)
            cDict    = cD.getScopDict()
            for k,v in cDict.items():
                self.__lfh.write("SCOP ACCESSION |%s|\n" % k)
                for k1,v1 in v.items():
                    self.__lfh.write("       %s value=%s\n" % (k1,v1))
            #
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()            

    def testGetScopDictionaryMethods(self): 
        """ Test dictionary methods -- 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            cD=ScopDictionary(configLevel='DEV',verbose=self.__verbose,log=self.__lfh)
            cDict    = cD.getScopDict()
            for k,v in cDict.items():
                self.__lfh.write("SCOP ACCESSION |%s|\n" % k)
                self.__lfh.write("Class:        %s|%s|%s\n" % (cD.getClassDescription(k)))
                self.__lfh.write("Fold:         %s|%s|%s\n" % (cD.getFoldDescription(k)))
                self.__lfh.write("Superfamily:  %s|%s|%s\n" % (cD.getSuperFamilyDescription(k)))
                self.__lfh.write("Family:       %s|%s|%s\n" % (cD.getFamilyDescription(k)))
                self.__lfh.write("Domain:       %s|%s|%s\n" % (cD.getDomainDescription(k)))                                
                
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()            



            
def suite():
    return unittest.makeSuite(ScopDictionaryTests,'test')


if __name__ == '__main__':
    unittest.main()
    
