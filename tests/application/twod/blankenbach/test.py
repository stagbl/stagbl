from stagbl_test_utils import create_stagbl_test
def test() :
    t = create_stagbl_test(__file__,
            'stagbldemo2d',
            1,
            [
                '--mode','blankenbach'
            ]
            )

    def comparefunc(t) :
        pass
        # TODO

    t.appendKeywords('%')
    t.appendKeywords('[')
    t.appendKeywords(']')
    t.appendKeywords('x')
    t.setVerifyMethod(comparefunc)

    return(t)
