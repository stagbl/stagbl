from stagbl_test_utils import create_stagbl_test
def test() :
    t = create_stagbl_test(__file__,
            'stagbldemo2d',
            9,
            [
                '-debug_ascii_dump',
                '-pc_type lu','-pc_factor_mat_solver_type mumps','-mat_mumps_icntl_23 100',
            ],
            'x.matlabascii.txt'
            )

    def comparefunc(t) :
        """ compare all lines not skipped by 'keywords' below """
        t.compareFloatingPointRelative('',1e-8,1e-11)

    t.appendKeywords('%')
    t.appendKeywords('[')
    t.appendKeywords(']')
    t.appendKeywords('x')
    t.setVerifyMethod(comparefunc)

    return(t)
