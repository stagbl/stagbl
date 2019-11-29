import os
from sciath.test import Test

def test() :
    # Note that we create this test in a special way, not with the usual
    # use of create_stagbl_test(). Thus, do not use this as an example test file.

    this_dir = os.path.split(os.path.abspath(__file__))[0]

    STAGBL_DIR = os.getenv('STAGBL_DIR')
    if not STAGBL_DIR :
        raise RuntimeError('STAGBL_DIR not defined in environment')

    STAGBL_ARCH = os.getenv('STAGBL_ARCH')
    if not STAGBL_ARCH :
        raise RuntimeError('STAGBL_ARCH not defined in environment')

    command_details = "stagbldemo2d"

    # Search the README.md to make sure the expected command is there
    with open(os.path.join(STAGBL_DIR,'README.rst')) as readme_file :
        if not command_details in readme_file.read():
            raise Exception("Did not find the expected quickstart command in README.rst! : " + command)

    command_full = os.path.join(STAGBL_DIR,STAGBL_ARCH,'bin',command_details)

    t = Test("quickstart1",1,command_full,os.path.join(this_dir,'expected'))
    t.setUseSandbox()

    def comparefunc(t) :
        # Check nothing
        pass

    t.setVerifyMethod(comparefunc)

    return(t)
