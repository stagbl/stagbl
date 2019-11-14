from stagbl_test_utils import create_stagbl_test

def test() :
    t = create_stagbl_test(__file__,
            'test_dmstag_preallocate',
            1,
            ['-preallocate 1','-log_view','-pc_type','lu','-pc_factor_mat_solver_type','umfpack'],
            )

    # Check that time doesn't change too much
    def comparefunc(t) :
        tol = 1.0e-1
        t.compareFloatingPointAbsolute('1:        Creation:',tol)

    t.setVerifyMethod(comparefunc)

    return(t)
