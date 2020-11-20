from .mapping import Mapping

class IdentityMapping(Mapping):
    """
    Represents an identity 1D/2D/3D Mapping object.

    Examples

    """
    _expressions = {'x': 'x1',
                    'y': 'x2',
                    'z': 'x3'}

#==============================================================================
class AffineMapping(Mapping):
    """
    Represents a 1D/2D/3D Affine Mapping object.

    Examples

    """
    _expressions = {'x': 'c1 + a11*x1 + a12*x2 + a13*x3',
                    'y': 'c2 + a21*x1 + a22*x2 + a23*x3',
                    'z': 'c3 + a31*x1 + a32*x2 + a33*x3'}

#==============================================================================
class PolarMapping(Mapping):
    """
    Represents a Polar 2D Mapping object (Annulus).

    Examples

    """
    _expressions = {'x': 'c1 + (rmin*(1-x1)+rmax*x1)*cos(x2)',
                    'y': 'c2 + (rmin*(1-x1)+rmax*x1)*sin(x2)'}

    _ldim        = 2
    _pdim        = 2

#==============================================================================
class TargetMapping(Mapping):
    """
    Represents a Target 2D Mapping object.

    Examples

    """
    _expressions = {'x': 'c1 + (1-k)*x1*cos(x2) - D*x1**2',
                    'y': 'c2 + (1+k)*x1*sin(x2)'}

    _ldim        = 2
    _pdim        = 2

#==============================================================================
class CzarnyMapping(Mapping):
    """
    Represents a Czarny 2D Mapping object.

    Examples

    """
    _expressions = {'x': '(1 - sqrt( 1 + eps*(eps + 2*x1*cos(x2)) )) / eps',
                    'y': 'c2 + (b / sqrt(1-eps**2/4) * x1 * sin(x2)) /'
                        '(2 - sqrt( 1 + eps*(eps + 2*x1*cos(x2)) ))'}

    _ldim        = 2
    _pdim        = 2

#==============================================================================
class CollelaMapping2D(Mapping):
    """
    Represents a Collela 2D Mapping object.

    Examples

    """
    _expressions = {'x': '2.*(x1 + eps*sin(2.*pi*k1*x1)*sin(2.*pi*k2*x2)) - 1.',
                    'y': '2.*(x2 + eps*sin(2.*pi*k1*x1)*sin(2.*pi*k2*x2)) - 1.'}

    _ldim        = 2
    _pdim        = 2

#==============================================================================
class TorusMapping(Mapping):
    """
    Represents a Torus 3D Mapping object.

    Examples

    """
    _expressions = {'x': '(R0+x1*cos(x2))*cos(x3)',
                    'y': '(R0+x1*cos(x2))*sin(x3)',
                    'z': 'x1*sin(x2)'}

    _ldim        = 2
    _pdim        = 2

#==============================================================================
class TwistedTargetMapping(Mapping):
    """
    Represents a Twisted Target 3D Mapping object.

    Examples

    """
    _expressions = {'x': 'c1 + (1-k)*x1*cos(x2) - D*x1**2',
                    'y': 'c2 + (1+k)*x1*sin(x2)',
                    'z': 'c3 + x3*x1**2*sin(2*x2)'}

    _ldim        = 2
    _pdim        = 2

