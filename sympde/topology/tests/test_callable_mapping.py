import numpy as np

from sympde.topology.mapping            import Mapping, BasicCallableMapping
from sympde.topology.analytical_mapping import IdentityMapping, AffineMapping
from sympde.topology.analytical_mapping import PolarMapping

# Tolerance for testing float equality
RTOL = 1e-15
ATOL = 1e-15

#==============================================================================
def test_identity_mapping_1d():

    F = IdentityMapping('F', dim=1)
    f = F.get_callable_mapping()

    assert f.ldim == 1
    assert f.pdim == 1

    x1_pts = [-0.7, 0.5, 33]

    for x1 in x1_pts:
        assert f(x1)              == (x1,)
        assert f.jacobian(x1)     == 1.0
        assert f.jacobian_inv(x1) == 1.0
        assert f.metric(x1)       == 1.0
        assert f.metric_det(x1)   == 1.0


def test_identity_mapping_2d():

    F = IdentityMapping('F', dim=2)
    f = F.get_callable_mapping()

    assert f.ldim == 2
    assert f.pdim == 2

    x1_pts = [-0.7, 0.5, 33]
    x2_pts = [-0.2, 1.3, 14]

    I = [[1, 0],
         [0, 1]]

    for x1 in x1_pts:
        for x2 in x2_pts:
            assert f(x1, x2) == (x1, x2)
            assert np.array_equal(f.jacobian    (x1, x2), I)
            assert np.array_equal(f.jacobian_inv(x1, x2), I)
            assert np.array_equal(f.metric      (x1, x2), I)
            assert f.metric_det(x1, x2) == 1.0


def test_identity_mapping_3d():

    F = IdentityMapping('F', dim=3)
    f = F.get_callable_mapping()

    assert f.ldim == 3
    assert f.pdim == 3

    x1_pts = [-0.5, 3]
    x2_pts = [-0.2, 2]
    x3_pts = [-1, 4.8]

    I = [[1, 0, 0],
         [0, 1, 0],
         [0, 0, 1]]

    for x1 in x1_pts:
        for x2 in x2_pts:
            for x3 in x3_pts:
                assert f(x1, x2, x3) == (x1, x2, x3)
                assert np.array_equal(f.jacobian    (x1, x2, x3), I)
                assert np.array_equal(f.jacobian_inv(x1, x2, x3), I)
                assert np.array_equal(f.metric      (x1, x2, x3), I)
                assert f.metric_det(x1, x2, x3) == 1.0


#------------------------------------------------------------------------------
def test_affine_mapping_1d():

    # x = 1 - x1
    params = {'c1': 1, 'a11': -1}

    F = AffineMapping('F', **params, dim=1)
    f = F.get_callable_mapping()

    assert f.ldim == 1
    assert f.pdim == 1

    assert f(0  ) == (1  ,)
    assert f(0.5) == (0.5,)
    assert f(1  ) == (0  ,)

    for x1 in [0, 0.5, 1]:
        assert f.jacobian(x1)     == -1
        assert f.jacobian_inv(x1) == -1
        assert f.metric(x1)       ==  1
        assert f.metric_det(x1)   ==  1


def test_affine_mapping_2d():

    c1, c2 = (1, 2)
    J = [[a11, a12], [a21, a22]] = [[3, 5], [-2, 4]]

    params = dict(c1=c1, c2=c2, a11=a11, a12=a12, a21=a21, a22=a22)

    F = AffineMapping('F', **params, dim=2)
    f = F.get_callable_mapping()

    assert f.ldim == 2
    assert f.pdim == 2

    J_det = a11 * a22 - a12 * a21
    J_inv = [[ a22 / J_det, -a12 / J_det],
             [-a21 / J_det,  a11 / J_det]]

    G = [[a11 * a11 + a21 * a21, a11 * a12 + a21 * a22],
         [a12 * a11 + a22 * a21, a12 * a12 + a22 * a22]]
    G_det = G[0][0] * G[1][1] - G[0][1] * G[1][0]

    x1_pts = [-0.7, 0.5, 33]
    x2_pts = [-0.2, 1.3, 14]

    for x1 in x1_pts:
        for x2 in x2_pts:
            x = c1 + a11 * x1 + a12 * x2
            y = c2 + a21 * x1 + a22 * x2
            assert f(x1, x2) == (x, y)
            assert np.array_equal(f.jacobian    (x1, x2), J    )
            assert np.array_equal(f.jacobian_inv(x1, x2), J_inv)
            assert np.array_equal(f.metric      (x1, x2), G    )
            assert f.metric_det(x1, x2) == G_det


def test_affine_mapping_3d():

    c1, c2, c3 = (-3, 1, 5)

    [a11, a12, a13] = [ 2,  7, -1]
    [a21, a22, a23] = [ 0,  3,  5]
    [a31, a32, a33] = [-4, -1,  1]

    params = dict(c1=c1, c2=c2, c3=c3,
            a11=a11, a12=a12, a13=a13,
            a21=a21, a22=a22, a23=a23,
            a31=a31, a32=a32, a33=a33,)

    F = AffineMapping('F', **params, dim=3)
    f = F.get_callable_mapping()

    assert f.ldim == 3
    assert f.pdim == 3

    J = np.array([[a11, a12, a13],
               [a21, a22, a23],
               [a31, a32, a33]])

    J_adj = np.array([[a22*a33-a23*a32, a13*a32-a12*a33, a12*a23-a13*a22],
                      [a23*a31-a21*a33, a11*a33-a13*a31, a13*a21-a11*a23],
                      [a21*a32-a22*a31, a12*a31-a11*a32, a11*a22-a12*a21]])

    J_det = (J * J_adj).sum()
    J_inv = J_adj / J_det

    G     = J.T @ J
    G_det = J_det**2

    x1_pts = [-0.5, 3]
    x2_pts = [-0.2, 2]
    x3_pts = [-1,   4]

    for x1 in x1_pts:
        for x2 in x2_pts:
            for x3 in x3_pts:
                x = c1 + a11 * x1 + a12 * x2 + a13 * x3
                y = c2 + a21 * x1 + a22 * x2 + a23 * x3
                z = c3 + a31 * x1 + a32 * x2 + a33 * x3
                assert np.allclose(f(x1, x2, x3), (x, y, z), rtol=RTOL, atol=ATOL)
                assert np.array_equal(f.jacobian    (x1, x2, x3), J    )
                assert np.array_equal(f.jacobian_inv(x1, x2, x3), J_inv)
                assert np.array_equal(f.metric      (x1, x2, x3), G    )
                assert f.metric_det(x1, x2, x3) == G_det


#------------------------------------------------------------------------------
def test_polar_mapping():

    c1, c2 = (-1, -1)
    rmin, rmax = (1, 2)

    params = dict(c1=c1, c2=c2, rmin=rmin, rmax=rmax)

    F = PolarMapping('F', **params)
    f = F.get_callable_mapping()

    assert f.ldim == 2
    assert f.pdim == 2

    assert f(0, 0) == (0, -1)
    assert f(1, 0) == (1, -1)
    assert np.allclose(f(0, np.pi/2), (-1, 0))

    x1_pts = [0, 0.5, 1]
    x2_pts = [-0.2, 1.3, 14]

    for x1 in x1_pts:
        for x2 in x2_pts:
            x  = c1 + (rmin * (1-x1) + rmax * x1) * np.cos(x2)
            y  = c2 + (rmin * (1-x1) + rmax * x1) * np.sin(x2)

            J  = [[(rmax-rmin) * np.cos(x2), -(rmin*(1-x1)+rmax*x1) * np.sin(x2)],
                  [(rmax-rmin) * np.sin(x2),  (rmin*(1-x1)+rmax*x1) * np.cos(x2)]]
            G  = [[(rmax-rmin)**2,                        0],
                  [             0, (rmin*(1-x1)+rmax*x1)**2]]

            G_det = G[0][0] * G[1][1] - G[0][1] * G[1][0]
            J_det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
            J_inv = [[ J[1][1] / J_det, -J[0][1] / J_det],
                     [-J[1][0] / J_det,  J[0][0] / J_det]]

            assert f(x1, x2) == (x, y)
            assert np.allclose(f.jacobian    (x1, x2), J    , rtol=RTOL, atol=ATOL)
            assert np.allclose(f.jacobian_inv(x1, x2), J_inv, rtol=RTOL, atol=ATOL)
            assert np.allclose(f.metric      (x1, x2), G    , rtol=RTOL, atol=ATOL)
            assert np.allclose(f.metric_det  (x1, x2), G_det, rtol=RTOL, atol=ATOL)


#------------------------------------------------------------------------------
# TODO [YG 23.02.2022]: add unit tests for other mappings

def test_identity_mapping_array_1d():

    F = IdentityMapping('F', dim=1)
    f = F.get_callable_mapping()

    assert f.ldim == 1
    assert f.pdim == 1

    x1 = np.array([-0.7, 0.5, 33])

    assert np.array_equal(f(x1)[0], x1)
    assert np.array_equal(f.jacobian(x1), np.ones((1, 1, 3)))
    assert np.array_equal(f.jacobian_inv(x1),np.ones((1, 1, 3)))
    assert np.array_equal(f.metric(x1), np.ones((1, 1, 3)))
    assert np.array_equal(f.metric_det(x1), np.ones(3))

def test_identity_mapping_array_2d():

    F = IdentityMapping('F', dim=2)
    f = F.get_callable_mapping()

    assert f.ldim == 2
    assert f.pdim == 2

    x1 = np.array([-0.7, 0.5, 33])
    x2 = np.array([-0.2, 1.3, 14])

    xx1, xx2 = np.meshgrid(x1, x2)

    I = np.zeros((2, 2, 3, 3)) + np.eye(2)[..., None, None]

    r1, r2 = f(xx1, xx2)

    assert np.array_equal(r1, xx1)
    assert np.array_equal(r2, xx2)
    assert np.array_equal(f.jacobian    (xx1, xx2), I)
    assert np.array_equal(f.jacobian_inv(xx1, xx2), I)
    assert np.array_equal(f.metric      (xx1, xx2), I)
    assert np.array_equal(f.metric_det  (xx1, xx2), np.ones((3, 3)))


def test_identity_mapping_array_3d():

    F = IdentityMapping('F', dim=3)
    f = F.get_callable_mapping()

    assert f.ldim == 3
    assert f.pdim == 3

    x1 = np.array([-0.5, 3])
    x2 = np.array([-0.2, 2])
    x3 = np.array([-1, 4.8])

    I = np.zeros((3, 3, 2, 2, 2)) + np.eye(3)[..., None, None, None]

    xx1, xx2, xx3 = np.meshgrid(x1, x2, x3)

    r1, r2, r3 = f(xx1, xx2, xx3)

    assert np.array_equal(r1, xx1)
    assert np.array_equal(r2, xx2)
    assert np.array_equal(r3, xx3)
    assert np.array_equal(f.jacobian    (xx1, xx2, xx3), I)
    assert np.array_equal(f.jacobian_inv(xx1, xx2, xx3), I)
    assert np.array_equal(f.metric      (xx1, xx2, xx3), I)
    assert np.array_equal(f.metric_det  (xx1, xx2, xx3), np.ones((2, 2, 2)))


#------------------------------------------------------------------------------
def test_affine_mapping_array_1d():

    # x = 1 - x1
    params = {'c1': 1, 'a11': -1}

    F = AffineMapping('F', **params, dim=1)
    f = F.get_callable_mapping()

    assert f.ldim == 1
    assert f.pdim == 1

    x1 = np.array([0, 0.5, 1])

    assert np.array_equal(f(x1)[0], [1, 0.5, 0])
    assert np.array_equal(f.jacobian(x1), [[[-1, -1, -1]]])
    assert np.array_equal(f.jacobian_inv(x1), [[[-1, -1, -1]]])
    assert np.array_equal(f.metric(x1), [[[1, 1, 1]]])
    assert np.array_equal(f.metric_det(x1), [1, 1, 1])


def test_affine_mapping_array_2d():

    c1, c2 = (1, 2)
    J = [[a11, a12], [a21, a22]] = [[3, 5], [-2, 4]]

    params = dict(c1=c1, c2=c2, a11=a11, a12=a12, a21=a21, a22=a22)

    F = AffineMapping('F', **params, dim=2)
    f = F.get_callable_mapping()

    assert f.ldim == 2
    assert f.pdim == 2

    J_det = a11 * a22 - a12 * a21
    J_inv = [[ a22 / J_det, -a12 / J_det],
             [-a21 / J_det,  a11 / J_det]]

    G = [[a11 * a11 + a21 * a21, a11 * a12 + a21 * a22],
         [a12 * a11 + a22 * a21, a12 * a12 + a22 * a22]]
    G_det = G[0][0] * G[1][1] - G[0][1] * G[1][0]

    x1 = np.array([-0.7, 0.5, 33])
    x2 = np.array([-0.2, 1.3, 14])

    xx1, xx2 = np.meshgrid(x1, x2)

    x = c1 + a11 * xx1 + a12 * xx2
    y = c2 + a21 * xx1 + a22 * xx2

    r1, r2 = f(xx1, xx2)

    J = np.zeros((2, 2, 3, 3)) + np.array(J)[..., None, None]
    J_inv = np.zeros((2, 2, 3, 3)) + np.array(J_inv)[..., None, None]
    G = np.zeros((2, 2, 3, 3)) + np.array(G)[..., None, None]
    G_det = np.full((3, 3), G_det)

    assert np.array_equal(r1, x)
    assert np.array_equal(r2, y)
    assert np.array_equal(f.jacobian    (xx1, xx2), J    )
    assert np.array_equal(f.jacobian_inv(xx1, xx2), J_inv)
    assert np.array_equal(f.metric      (xx1, xx2), G    )
    assert np.array_equal(f.metric_det  (xx1, xx2), G_det)


def test_affine_mapping_array_3d():

    c1, c2, c3 = (-3, 1, 5)

    [a11, a12, a13] = [ 2,  7, -1]
    [a21, a22, a23] = [ 0,  3,  5]
    [a31, a32, a33] = [-4, -1,  1]

    params = dict(c1=c1, c2=c2, c3=c3,
            a11=a11, a12=a12, a13=a13,
            a21=a21, a22=a22, a23=a23,
            a31=a31, a32=a32, a33=a33,)

    F = AffineMapping('F', **params, dim=3)
    f = F.get_callable_mapping()

    assert f.ldim == 3
    assert f.pdim == 3

    J = np.array([[a11, a12, a13],
                  [a21, a22, a23],
                  [a31, a32, a33]])

    J_adj = np.array([[a22*a33-a23*a32, a13*a32-a12*a33, a12*a23-a13*a22],
                      [a23*a31-a21*a33, a11*a33-a13*a31, a13*a21-a11*a23],
                      [a21*a32-a22*a31, a12*a31-a11*a32, a11*a22-a12*a21]])

    J_det = (J * J_adj).sum()
    J_inv = J_adj / J_det

    G     = J.T @ J
    G_det = J_det**2

    x1 = np.array([-0.5, 3])
    x2 = np.array([-0.2, 2])
    x3 = np.array([-1,   4])

    xx1, xx2, xx3 = np.meshgrid(x1, x2, x3)

    r1, r2, r3 = f(xx1, xx2, xx3)

    x = c1 + a11 * xx1 + a12 * xx2 + a13 * xx3
    y = c2 + a21 * xx1 + a22 * xx2 + a23 * xx3
    z = c3 + a31 * xx1 + a32 * xx2 + a33 * xx3

    J = np.zeros((3, 3, 2, 2, 2)) + J[..., None, None, None]
    J_inv = np.zeros((3, 3, 2, 2, 2)) + J_inv[..., None, None, None]
    G = np.zeros((3, 3, 2, 2, 2)) + G[..., None, None, None]
    G_det = np.full((2, 2, 2), G_det)

    assert np.allclose(r1, x, rtol=RTOL, atol=ATOL)
    assert np.allclose(r2, y, rtol=RTOL, atol=ATOL)
    assert np.allclose(r3, z, rtol=RTOL, atol=ATOL)
    assert np.array_equal(f.jacobian    (xx1, xx2, xx3), J    )
    assert np.array_equal(f.jacobian_inv(xx1, xx2, xx3), J_inv)
    assert np.array_equal(f.metric      (xx1, xx2, xx3), G    )
    assert np.array_equal(f.metric_det  (xx1, xx2, xx3), G_det)


#------------------------------------------------------------------------------
def test_polar_mapping_array():

    c1, c2 = (-1, -1)
    rmin, rmax = (1, 2)

    params = dict(c1=c1, c2=c2, rmin=rmin, rmax=rmax)

    F = PolarMapping('F', **params)
    f = F.get_callable_mapping()

    assert f.ldim == 2
    assert f.pdim == 2

    assert f(0, 0) == (0, -1)
    assert f(1, 0) == (1, -1)
    assert np.allclose(f(0, np.pi/2), (-1, 0))

    x1 = np.array([0, 0.5, 1])
    x2 = np.array([-0.2, 1.3, 14])

    xx1, xx2 = np.meshgrid(x1, x2)

    r1, r2 = f(xx1, xx2)

    x  = c1 + (rmin * (1-xx1) + rmax * xx1) * np.cos(xx2)
    y  = c2 + (rmin * (1-xx1) + rmax * xx1) * np.sin(xx2)

    J  = np.array([[(rmax-rmin) * np.cos(xx2), -(rmin*(1-xx1)+rmax*xx1) * np.sin(xx2)],
                   [(rmax-rmin) * np.sin(xx2),  (rmin*(1-xx1)+rmax*xx1) * np.cos(xx2)]])

    G  = np.array([[np.full_like(xx1, (rmax-rmin)**2),         np.zeros_like(xx1)],
                   [               np.zeros_like(xx1), (rmin*(1-xx1)+rmax*xx1)**2]])

    G_det = G[0, 0] * G[1, 1] - G[0, 1] * G[1, 0]
    J_det = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    J_inv = np.array([[ J[1, 1] / J_det, -J[0, 1] / J_det],
                      [-J[1, 0] / J_det,  J[0, 0] / J_det]])

    assert np.array_equal(r1, x)
    assert np.array_equal(r2, y)
    assert np.allclose(f.jacobian    (xx1, xx2), J    , rtol=RTOL, atol=ATOL)
    assert np.allclose(f.jacobian_inv(xx1, xx2), J_inv, rtol=RTOL, atol=ATOL)
    assert np.allclose(f.metric      (xx1, xx2), G    , rtol=RTOL, atol=ATOL)
    assert np.allclose(f.metric_det  (xx1, xx2), G_det, rtol=RTOL, atol=ATOL)

#==============================================================================
def test_user_defined_callable_mapping():

    class UserIdentity(BasicCallableMapping):
        """ Identity in N dimensions.
        """

        def __init__(self, ndim):
            self._ndim = ndim

        def __call__(self, *eta):
            assert len(eta) == self._ndim
            return eta

        def jacobian(self, *eta):
            assert len(eta) == self._ndim
            return np.eye(self._ndim)

        def jacobian_inv(self, *eta):
            assert len(eta) == self._ndim
            return np.eye(self._ndim)

        def metric(self, *eta):
            assert len(eta) == self._ndim
            return np.eye(self._ndim)

        def metric_det(self, *eta):
            assert len(eta) == self._ndim
            return 1.0

        @property
        def ldim(self):
            return self._ndim

        @property
        def pdim(self):
            return self._ndim

    F = Mapping('F', ldim = 3, pdim = 3) # Declare undefined symbolic mapping
    f = UserIdentity(3)        # Create user-defined callable mapping
    F.set_callable_mapping(f)  # Attach callable mapping to symbolic mapping

    assert F.get_callable_mapping() is f
    assert f(4, 5, 6) == (4, 5, 6)
    assert np.array_equal(f.jacobian(1, 2, 3)     , [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert np.array_equal(f.jacobian_inv(-7, 0, 7), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert np.array_equal(f.metric(0, 1, 1)       , [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert f.metric_det(-5, 8, -9) == 1.0
    assert f.ldim == 3
    assert f.pdim == 3
