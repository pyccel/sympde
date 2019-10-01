from sympde.core.utils import expand_name_patterns

#==============================================================================
def test_expand_name_patterns():

    assert expand_name_patterns('x,y,z') == ('x', 'y', 'z')
    assert expand_name_patterns('x y z') == ('x', 'y', 'z')

    assert expand_name_patterns('x,') == ('x',)
    assert expand_name_patterns('x', seq=True) == ('x',)

    assert expand_name_patterns(('a', 'b', 'c')) == ('a', 'b', 'c')
    assert expand_name_patterns(['a', 'b', 'c']) == ['a', 'b', 'c']
    assert expand_name_patterns({'a', 'b', 'c'}) == {'a', 'b', 'c'}

    assert expand_name_patterns('x:4') == ('x0', 'x1', 'x2', 'x3')
    assert expand_name_patterns('x5:10') == ('x5', 'x6', 'x7', 'x8', 'x9')
    assert expand_name_patterns('x7(1:3)') == ('x71', 'x72')
    assert expand_name_patterns('x2:5, y:2') == ('x2', 'x3', 'x4', 'y0', 'y1')
    assert expand_name_patterns(('x2:5', 'y:2')) == (('x2', 'x3', 'x4'), ('y0', 'y1'))

    assert expand_name_patterns(':c') == ('a', 'b', 'c')
    assert expand_name_patterns('x:c') == ()
    assert expand_name_patterns('x:z') == ('x', 'y', 'z')
    assert expand_name_patterns('x(:c)') == ('xa', 'xb', 'xc')
    assert expand_name_patterns('a:b, x:z') == ('a', 'b', 'x', 'y', 'z')
    assert expand_name_patterns(('a:b', 'x:z')) == (('a', 'b'), ('x', 'y', 'z'))

    assert expand_name_patterns('x:2(1:3)') == ('x01', 'x02', 'x11', 'x12')
    assert expand_name_patterns(':3:2') == ('00', '01', '10', '11', '20', '21')

    assert expand_name_patterns('x((a:b))') == ('x(a)', 'x(b)')
    assert expand_name_patterns(r'x(:1\,:2)') == ('x(0,0)', 'x(0,1)')
