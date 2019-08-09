from sympde.topology import InteriorDomain, Union
from sympde.topology import Boundary, NormalVector, TangentVector
from sympde.topology import Connectivity, Edge
from sympde.topology import Domain, ElementDomain
from sympde.topology import Area
from sympde.topology import Domain, Line, Square, Cube

def two_patches_Square():

    A = Square('A')
    B = Square('B')

    A = A.interior
    B = B.interior
    
    connectivity = Connectivity()

    bnd_A_1 = Boundary('Gamma_1', A, axis=0, ext=-1)
    bnd_A_2 = Boundary('Gamma_2', A, axis=0, ext=1)
    bnd_A_3 = Boundary('Gamma_3', A, axis=1, ext=-1)
    bnd_A_4 = Boundary('Gamma_4', A, axis=1, ext=1)

    bnd_B_1 = Boundary('Gamma_1', B, axis=0, ext=-1)
    bnd_B_2 = Boundary('Gamma_2', B, axis=0, ext=1)
    bnd_B_3 = Boundary('Gamma_3', B, axis=1, ext=-1)
    bnd_B_4 = Boundary('Gamma_4', B, axis=1, ext=1)

    connectivity['I'] = (bnd_A_2, bnd_B_1)

    Omega = Domain('Omega',
                   interiors=[A, B],
                   boundaries=[bnd_A_1, bnd_A_2, bnd_A_3, bnd_A_4, bnd_B_1, bnd_B_2, bnd_B_3, bnd_B_4],
                   connectivity=connectivity)

    return Omega
