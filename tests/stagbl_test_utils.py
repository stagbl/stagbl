import os
import pyTestHarness.test

def create_stagbl_test(calling_file,exec_name,ranks,arguments) :
    """ Set up a StagBL test based on a the absolute path of the calling file,
        the name of an executable in StagBL's bin directory,
        and a list of arguments """

    # Require STAGBL_DIR and STAGBL_ARCH in the environment
    STAGBL_DIR = os.getenv('STAGBL_DIR')
    if not STAGBL_DIR :
        raise RuntimeError('STAGBL_DIR not defined in environment')

    STAGBL_ARCH = os.getenv('STAGBL_ARCH')
    if not STAGBL_ARCH :
        raise RuntimeError('STAGBL_ARCH not defined in environment')

    # Construct the test name from the path to the directory defining the test
    calling_dir = os.path.split(os.path.abspath(calling_file))[0]
    test_name = os.path.relpath(calling_dir,os.path.join(STAGBL_DIR,'tests')).replace(os.sep,'.')

    # Use a default expected file
    expectedFile = os.path.join(calling_dir,'expected')

    # Assume a standard location for the executable and add space-separated arguments
    launch = os.path.join(STAGBL_DIR,STAGBL_ARCH,'bin',exec_name) + ' ' + ' '.join(arguments)

    # Create and set up a pth test, which can be operated on later
    t = pyTestHarness.test.Test(test_name,ranks,launch,expectedFile)
    t.setUseSandbox()

    return(t)
