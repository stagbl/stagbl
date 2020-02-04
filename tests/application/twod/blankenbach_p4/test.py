from stagbl_test_utils import create_stagbl_test

def test() :
    t = create_stagbl_test(__file__,
            'stagbldemo2d',
            4,
            [
                '-mode','blankenbach',
                '-nsteps','7',
                '-dt','1e17',
                '-stag_grid_x','39',
                '-stag_grid_y','37',
            ]
            )

    def comparefunc(t) :
        key = 'Nusselt number:'
        t.compareFloatingPointRelative(key,1.0e-5)
        key = 'Rayleigh number:'
        t.compareFloatingPointRelative(key,1.0e-5)

    t.setVerifyMethod(comparefunc)

    return(t)
