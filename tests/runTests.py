#!/usr/bin/env python
#################################################################################
#                                                                               #
#                                StagBL Tests                                   #
#                                                                               #
#################################################################################
# Usage:                                                                        #
#  ./runTests.py [-h] [-other_pyTestHarness_options]                            #
#                                                                               #
#################################################################################

import os
import sys
import json

thisDir = os.path.dirname(os.path.abspath(__file__))

# STAGBL_DIR and STAGBL_ARCH are used to make easily-copyable tests
STAGBL_DIR = os.getenv('STAGBL_DIR')
if not STAGBL_DIR :
    print('You must define STAGBL_DIR in your environment. Exiting.')
    sys.exit(1)
STAGBL_ARCH = os.getenv('STAGBL_ARCH')
if not STAGBL_ARCH :
    print('You must define STAGBL_ARCH in your environment. Exiting.')
    sys.exit(1)
if os.path.join(STAGBL_DIR,'tests') != thisDir :
    print('STAGBL_DIR is not set correctly in your environment. Exiting.')
    sys.exit(1)

# TODO load JSON

# TODO update "bin" to something from that json 

chkpath = os.path.join(STAGBL_DIR,STAGBL_ARCH,'bin')
if not os.path.exists(chkpath) :
    print(chkpath,' not found. Set STAGBL_DIR and STAGBL_ARCH properly. Exiting.')
    sys.exit(1)

# bitbucket.org/dmay/pythontestharness
sys.path.append(os.path.join(STAGBL_DIR,'pythontestharness','lib'))  # overrides
try :
  import pyTestHarness.harness as pth_harness
  import pyTestHarness.test as pth_test
except Exception as errorMessage :
  if not sys.exc_info()[-1].tb_next :     # Check that the traceback has depth 1
    traceback.print_exc()
    print('********************')
    print('The required python library pyTestHarness was not found. Exiting.')
    print('If pyTestHarness is installed on your system, ensure pythontestharness/lib is included in the environment variable PYTHONPATH.')
    print('If pyTestHarness is not installed, obtain the source by executing the following:')
    print('  git clone https://bitbucket.org/dmay/pythontestharness ' + os.path.join(STAGBL_DIR,'pythontestharness'))
    print('********************')
    sys.exit(1)
  raise

# Define a subclass of pyTestHarness.Test
class StagBLTest(pth_test.Test) :
    def __init__(self) :
        pth_test.Test.__init__(self)

#  Interpret any child directory containing 'test.py' as defining as a test
#  with a name defined as the relative path to the containing directory.
allTests = []
for (root, dirs, files) in os.walk(thisDir) :
    if 'test.py' in files :
        relpath = os.path.relpath(root,thisDir)
        mod = relpath.replace(os.sep,'.') + '.test'
        exec('import ' + mod)
        exec('allTests.append('+mod+'.test())')

# Run tests
h = pth_harness.pthHarness(allTests)
h.execute()
h.verify()
