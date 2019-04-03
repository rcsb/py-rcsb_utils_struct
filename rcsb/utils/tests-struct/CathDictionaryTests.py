##
# File:    CathDictionaryTests.py
# Author:  jdw
# Date:    9-May-2011
# Version: 0.001
#
# Updates:
#
##
"""
Test cases for operations that read CATH term and class data from flat files -

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.001"

import sys, unittest, traceback
import sys, time, os, os.path, shutil


from sbkb.utils.ConfigInfo             import ConfigInfo
from sbkb.readers.cath.CathDictionary  import CathDictionary

class CathDictionaryTests(unittest.TestCase):
    def setUp(self):
        self.__verbose=True
        self.__lfh=sys.stderr

    def tearDown(self):
        pass

    def testGetCathTermDictionary(self): 
        """ Get dictionary of CATH terms 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            cD=CathDictionary(configLevel='DEV',verbose=self.__verbose,log=self.__lfh)
            cDict    = cD.getCathTermDict()
            for k,v in cDict.items():
                self.__lfh.write("CATH ACCESSION |%s|\n" % k)
                self.__lfh.write("CATH CODE:        %s\n" % cD.getCathCode(k))
                self.__lfh.write("CATH DESCRIPTION: %s\n" % cD.getCathText(k))
                #
                for k1,v1 in v.items():
                    self.__lfh.write("       %s value=%s\n" % (k1,v1))
            #
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()            

    def testGetCathDictionaryMethods(self): 
        """ Test dictionary methods -- 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            cD=CathDictionary(configLevel='DEV',verbose=self.__verbose,log=self.__lfh)
            cDict    = cD.getCathTermDict()
            for k,v in cDict.items():
                self.__lfh.write("%s\t" % str(k).lower())
                self.__lfh.write("%s\t" % cD.getCathCode(k))        
                self.__lfh.write("%s\n" % cD.getCathTextTdd(k))
                self.__lfh.write("%s|%s\n" % cD.getCathHomology(k))
                self.__lfh.write("%s|%s\n" % cD.getCathTopology(k))
                self.__lfh.write("%s|%s\n" % cD.getCathArchitecture(k))
                self.__lfh.write("%s|%s\n" % cD.getCathClass(k))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()            



            
def suite():
    return unittest.makeSuite(CathDictionaryTests,'test')


if __name__ == '__main__':
    unittest.main()
    
