from .base_analytic_mapping import BaseAnalyticMapping

class IdentityMapping(BaseAnalyticMapping):
    """
    Represents an identity 1D/2D/3D BaseAnalyticMapping object.

    Examples

    """
    _expressions = {'x': 'x1',
                    'y': 'x2',
                    'z': 'x3'}

#==============================================================================
class AffineMapping(BaseAnalyticMapping):
    """
    Represents a 1D/2D/3D Affine BaseAnalyticMapping object.

    Examples

    """
    _expressions = {'x': 'c1 + a11*x1 + a12*x2 + a13*x3',
                    'y': 'c2 + a21*x1 + a22*x2 + a23*x3',
                    'z': 'c3 + a31*x1 + a32*x2 + a33*x3'}

#==============================================================================
class PolarMapping(BaseAnalyticMapping):
    """
    Represents a Polar 2D BaseAnalyticMapping object (Annulus).

    Examples

    """
    _expressions = {'x': 'c1 + (rmin*(1-x1)+rmax*x1)*cos(x2)',
                    'y': 'c2 + (rmin*(1-x1)+rmax*x1)*sin(x2)'}

    _ldim        = 2
    _pdim        = 2

#==============================================================================
class TargetMapping(BaseAnalyticMapping):
    """
    Represents a Target 2D BaseAnalyticMapping object.

    Examples

    """
    _expressions = {'x': 'c1 + (1-k)*x1*cos(x2) - D*x1**2',
                    'y': 'c2 + (1+k)*x1*sin(x2)'}

    _ldim        = 2
    _pdim        = 2

#==============================================================================
class CzarnyMapping(BaseAnalyticMapping):
    """
    Represents a Czarny 2D BaseAnalyticMapping object.

    Examples

    """
    _expressions = {'x': '(1 - sqrt( 1 + eps*(eps + 2*x1*cos(x2)) )) / eps',
                    'y': 'c2 + (b / sqrt(1-eps**2/4) * x1 * sin(x2)) /'
                        '(2 - sqrt( 1 + eps*(eps + 2*x1*cos(x2)) ))'}

    _ldim        = 2
    _pdim        = 2

#==============================================================================
class CollelaMapping2D(BaseAnalyticMapping):

    _expressions = {'x': '2.*(x1 + eps*sin(2.*pi*k1*x1)*sin(2.*pi*k2*x2)) - 1.',
                    'y': '2.*(x2 + eps*sin(2.*pi*k1*x1)*sin(2.*pi*k2*x2)) - 1.'}

    _ldim = 2
    _pdim = 2


#==============================================================================
class CollelaMapping3D(BaseAnalyticMapping):

    _expressions = {'x': 'k1*(x1 + eps*sin(2.*pi*x1)*sin(2.*pi*x2))',
                    'y': 'k2*(x2 + eps*sin(2.*pi*x1)*sin(2.*pi*x2))',
                    'z': 'k3*x3'}

    _ldim        = 3
    _pdim        = 3

#==============================================================================
class TorusMapping(BaseAnalyticMapping):
    """
    Parametrization of a torus (or a portion of it) of major radius R0, using
    toroidal coordinates (x1, x2, x3) = (r, theta, phi), where:

      - minor radius    0 <= r < R0
      - poloidal angle  0 <= theta < 2 pi
      - toroidal angle  0 <= phi < 2 pi

    """
    _expressions = {'x': '(R0 + x1 * cos(x2)) * cos(x3)',
                    'y': '(R0 + x1 * cos(x2)) * sin(x3)',
                    'z':       'x1 * sin(x2)'}

    _ldim        = 3
    _pdim        = 3

#==============================================================================
# TODO [YG, 07.10.2022]: add test in sympde/topology/tests/test_logical_expr.py
class TorusSurfaceMapping(BaseAnalyticMapping):
    """
    3D surface obtained by "slicing" the torus above at r = a.
    The parametrization uses the coordinates (x1, x2) = (theta, phi), where:

      - poloidal angle  0 <= theta < 2 pi
      - toroidal angle  0 <= phi < 2 pi

    """
    _expressions = {'x': '(R0 + a * cos(x1)) * cos(x2)',
                    'y': '(R0 + a * cos(x1)) * sin(x2)',
                    'z':       'a * sin(x1)'}

    _ldim        = 2
    _pdim        = 3

#==============================================================================
# TODO [YG, 07.10.2022]: add test in sympde/topology/tests/test_logical_expr.py
class TwistedTargetSurfaceMapping(BaseAnalyticMapping):
    """
    3D surface obtained by "twisting" the TargetMapping out of the (x, y) plane

    """
    _expressions = {'x': 'c1 + (1-k) * x1 * cos(x2) - D *x1**2',
                    'y': 'c2 + (1+k) * x1 * sin(x2)',
                    'z': 'c3 + x1**2 * sin(2*x2)'}

    _ldim        = 2
    _pdim        = 3

#==============================================================================
class TwistedTargetMapping(BaseAnalyticMapping):
    """
    3D volume obtained by "extruding" the TwistedTargetSurfaceMapping along z.

    """
    _expressions = {'x': 'c1 + (1-k) * x1 * cos(x2) - D * x1**2',
                    'y': 'c2 + (1+k) * x1 * sin(x2)',
                    'z': 'c3 + x3 * x1**2 * sin(2*x2)'}

    _ldim        = 3
    _pdim        = 3

#==============================================================================
class SphericalMapping(BaseAnalyticMapping):
    """
    Parametrization of a sphere (or a portion of it) using spherical
    coordinates (x1, x2, x3) = (r, theta, phi), where:

      - radius      r >= 0
      - inclination 0 <= theta <= pi
      - azimuth     0 <= phi < 2 pi

    """
    _expressions = {'x': 'x1 * sin(x2) * cos(x3)',
                    'y': 'x1 * sin(x2) * sin(x3)',
                    'z': 'x1 * cos(x2)'}

    _ldim        = 3
    _pdim        = 3

#==============================================================================
class Collela3D(BaseAnalyticMapping):

    _expressions = {'x':'2.*(x1 + 0.1*sin(2.*pi*x1)*sin(2.*pi*x2)) - 1.',
                    'y':'2.*(x2 + 0.1*sin(2.*pi*x1)*sin(2.*pi*x2)) - 1.',
                    'z':'2.*x3  - 1.'}

#==============================================================================   
class TransposedPolarMapping(BaseAnalyticMapping):
    """
    Represents a Transposed (x1 <> x2) Polar 2D Mapping object (Annulus).

    Examples

    """
    _expressions = {'x': 'c1 + (rmin*(1-x2)+rmax*x2)*cos(x1)',
                    'y': 'c2 + (rmin*(1-x2)+rmax*x2)*sin(x1)'}

    _ldim        = 2
    _pdim        = 2

#==============================================================================    
class TransposedPolarMapping(BaseAnalyticMapping):
    """
    Represents a Transposed (x1 <> x2) Polar 2D Mapping object (Annulus).

    Examples

    """
    _expressions = {'x': 'c1 + (rmin*(1-x2)+rmax*x2)*cos(x1)',
                    'y': 'c2 + (rmin*(1-x2)+rmax*x2)*sin(x1)'}

    _ldim        = 2
    _pdim        = 2
