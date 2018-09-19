import os
import pyTestHarness.unittest as pth

def test() :
    STAGBL_DIR = os.getenv('STAGBL_DIR')
    if not STAGBL_DIR :
        raise RuntimeError('STAGBL_DIR not defined in environment')
    STAGBL_ARCH = os.getenv('STAGBL_ARCH')
    if not STAGBL_ARCH :
        raise RuntimeError('STAGBL_ARCH not defined in environment')
    thisDir = os.path.split(os.path.abspath(__file__))[0]
    testName = os.path.relpath(thisDir,os.path.join(STAGBL_DIR,'tests')).replace(os.sep,'.')
    ranks = 2
    launch = os.path.join(STAGBL_DIR,STAGBL_ARCH,'bin','stagbldemo2d')
    expectedFile = os.path.join(thisDir,'expected')

    # Do nothing (test will almost certainly pass)
    def comparefunc(t) :
        pass

    t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
    t.setUseSandbox()
    t.setVerifyMethod(comparefunc)

    return(t)
