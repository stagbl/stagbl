from stagbl_test_utils import create_stagbl_test

def test() :
    t = create_stagbl_test(__file__,
            'stagbldemo3d',
            2,
            [
                '-stag_grid_x 5','-stag_grid_y 5','-stag_grid_z 5',
                '-pc_type lu','-pc_factor_mat_solver_package mumps',
            ]
            )

    # Do nothing (test will almost certainly pass)
    def comparefunc(t) :
        pass

    t.setVerifyMethod(comparefunc)

    return(t)
