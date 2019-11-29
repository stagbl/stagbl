from stagbl_test_utils import create_stagbl_test

def test() :
    t = create_stagbl_test(__file__,
            'stagbldemo2d',
            1,
            [],
            )

    # Do nothing (test will almost certainly pass)
    def comparefunc(t) :
        pass

    t.setVerifyMethod(comparefunc)

    return(t)
