from stagbl_test_utils import create_stagbl_test

def test() :
    t = create_stagbl_test(__file__,
            'stagbldemo3d',
            1,
            ['-stag_grid_x 5','-stag_grid_y 5','-stag_grid_z 5']
            )

    # Do nothing (test will almost certainly pass)
    def comparefunc(t) :
        pass

    t.setVerifyMethod(comparefunc)

    return(t)
