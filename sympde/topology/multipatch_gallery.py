# --------------------------------------------------------------------------- #
# This file is part of SymPDE.                                              #
# --------------------------------------------------------------------------- #
"""Gallery of symbolic 2D and 3D multipatch domains."""

from argparse import ArgumentParser
from pathlib import Path

import numpy as np

from .analytical_mapping import (
    AffineMapping,
    IdentityMapping,
    PolarMapping,
    TorusMapping,
    TransposedPolarMapping,
)
from .domain import Cube, Domain, Square

__all__ = (
    'build_annulus_3',
    'build_annulus_4',
    'build_cartesian_multipatch_domain_2d',
    'build_curved_l_shape',
    'build_multipatch_domain_2d',
    'build_multipatch_domain_3d',
    'build_oriented_two_cube',
    'build_pretzel',
    'build_pretzel_annulus',
    'build_pretzel_debug',
    'build_pretzel_f',
    'build_self_connected_patch',
    'build_square_2',
    'build_square_4',
    'build_square_6',
    'build_square_8',
    'build_square_9',
    'build_three_patch_shared_edge',
    'build_three_patch_shared_vertex',
    'build_torus_3d',
    'build_two_cube_3d',
    'build_two_patch_twisted_strip',
    'plot_multipatch_domain',
)


# =============================================================================
# 2D domains
# =============================================================================

def build_square_2():
    """Build a square decomposed into two horizontal patches."""
    logical_patch_1 = Square(
        'P1',
        bounds1=(0., np.pi),
        bounds2=(0., np.pi / 2),
    )
    logical_patch_2 = Square(
        'P2',
        bounds1=(0., np.pi),
        bounds2=(np.pi / 2, np.pi),
    )
    patch_1 = IdentityMapping('M1', dim=2)(logical_patch_1)
    patch_2 = IdentityMapping('M2', dim=2)(logical_patch_2)

    patches = (patch_1, patch_2)
    connectivity = (
        ((patch_1, 1, +1), (patch_2, 1, -1), +1),
    )
    return Domain.join(patches, connectivity, name='square_2')


def build_square_4():
    """Build a square decomposed into a two-by-two patch grid."""
    logical_patch_a = Square(
        'A',
        bounds1=(0., np.pi / 2),
        bounds2=(0., np.pi / 2),
    )
    logical_patch_b = Square(
        'B',
        bounds1=(np.pi / 2, np.pi),
        bounds2=(0., np.pi / 2),
    )
    logical_patch_c = Square(
        'C',
        bounds1=(0., np.pi / 2),
        bounds2=(np.pi / 2, np.pi),
    )
    logical_patch_d = Square(
        'D',
        bounds1=(np.pi / 2, np.pi),
        bounds2=(np.pi / 2, np.pi),
    )
    patch_a = IdentityMapping('M1', dim=2)(logical_patch_a)
    patch_b = IdentityMapping('M2', dim=2)(logical_patch_b)
    patch_c = IdentityMapping('M3', dim=2)(logical_patch_c)
    patch_d = IdentityMapping('M4', dim=2)(logical_patch_d)

    patches = (patch_a, patch_b, patch_c, patch_d)
    connectivity = (
        ((patch_a, 0, +1), (patch_b, 0, -1), +1),
        ((patch_a, 1, +1), (patch_c, 1, -1), +1),
        ((patch_c, 0, +1), (patch_d, 0, -1), +1),
        ((patch_b, 1, +1), (patch_d, 1, -1), +1),
    )
    return Domain.join(patches, connectivity, name='square_4')


def build_square_6():
    """Build a square decomposed into a two-by-three patch grid."""
    logical_patch_1 = Square(
        'P1',
        bounds1=(0., np.pi / 2),
        bounds2=(0., np.pi / 3),
    )
    logical_patch_2 = Square(
        'P2',
        bounds1=(np.pi / 2, np.pi),
        bounds2=(0., np.pi / 3),
    )
    logical_patch_3 = Square(
        'P3',
        bounds1=(0., np.pi / 2),
        bounds2=(np.pi / 3, 2 * np.pi / 3),
    )
    logical_patch_4 = Square(
        'P4',
        bounds1=(np.pi / 2, np.pi),
        bounds2=(np.pi / 3, 2 * np.pi / 3),
    )
    logical_patch_5 = Square(
        'P5',
        bounds1=(0., np.pi / 2),
        bounds2=(2 * np.pi / 3, np.pi),
    )
    logical_patch_6 = Square(
        'P6',
        bounds1=(np.pi / 2, np.pi),
        bounds2=(2 * np.pi / 3, np.pi),
    )
    patch_1 = IdentityMapping('M1', dim=2)(logical_patch_1)
    patch_2 = IdentityMapping('M2', dim=2)(logical_patch_2)
    patch_3 = IdentityMapping('M3', dim=2)(logical_patch_3)
    patch_4 = IdentityMapping('M4', dim=2)(logical_patch_4)
    patch_5 = IdentityMapping('M5', dim=2)(logical_patch_5)
    patch_6 = IdentityMapping('M6', dim=2)(logical_patch_6)

    patches = (patch_1, patch_2, patch_3, patch_4, patch_5, patch_6)
    connectivity = (
        ((patch_1, 0, +1), (patch_2, 0, -1), +1),
        ((patch_3, 0, +1), (patch_4, 0, -1), +1),
        ((patch_5, 0, +1), (patch_6, 0, -1), +1),
        ((patch_1, 1, +1), (patch_3, 1, -1), +1),
        ((patch_3, 1, +1), (patch_5, 1, -1), +1),
        ((patch_2, 1, +1), (patch_4, 1, -1), +1),
        ((patch_4, 1, +1), (patch_6, 1, -1), +1),
    )
    return Domain.join(patches, connectivity, name='square_6')


def build_square_8():
    """Build a three-by-three square grid without its center patch."""
    logical_patch_1 = Square(
        'P1',
        bounds1=(0., np.pi / 3),
        bounds2=(0., np.pi / 3),
    )
    logical_patch_2 = Square(
        'P2',
        bounds1=(np.pi / 3, 2 * np.pi / 3),
        bounds2=(0., np.pi / 3),
    )
    logical_patch_3 = Square(
        'P3',
        bounds1=(2 * np.pi / 3, np.pi),
        bounds2=(0., np.pi / 3),
    )
    logical_patch_4 = Square(
        'P4',
        bounds1=(0., np.pi / 3),
        bounds2=(np.pi / 3, 2 * np.pi / 3),
    )
    logical_patch_5 = Square(
        'P5',
        bounds1=(2 * np.pi / 3, np.pi),
        bounds2=(np.pi / 3, 2 * np.pi / 3),
    )
    logical_patch_6 = Square(
        'P6',
        bounds1=(0., np.pi / 3),
        bounds2=(2 * np.pi / 3, np.pi),
    )
    logical_patch_7 = Square(
        'P7',
        bounds1=(np.pi / 3, 2 * np.pi / 3),
        bounds2=(2 * np.pi / 3, np.pi),
    )
    logical_patch_8 = Square(
        'P8',
        bounds1=(2 * np.pi / 3, np.pi),
        bounds2=(2 * np.pi / 3, np.pi),
    )
    patch_1 = IdentityMapping('M1', dim=2)(logical_patch_1)
    patch_2 = IdentityMapping('M2', dim=2)(logical_patch_2)
    patch_3 = IdentityMapping('M3', dim=2)(logical_patch_3)
    patch_4 = IdentityMapping('M4', dim=2)(logical_patch_4)
    patch_5 = IdentityMapping('M5', dim=2)(logical_patch_5)
    patch_6 = IdentityMapping('M6', dim=2)(logical_patch_6)
    patch_7 = IdentityMapping('M7', dim=2)(logical_patch_7)
    patch_8 = IdentityMapping('M8', dim=2)(logical_patch_8)

    patches = (
        patch_1, patch_2, patch_3, patch_4,
        patch_5, patch_6, patch_7, patch_8,
    )
    connectivity = (
        ((patch_1, 0, +1), (patch_2, 0, -1), +1),
        ((patch_2, 0, +1), (patch_3, 0, -1), +1),
        ((patch_6, 0, +1), (patch_7, 0, -1), +1),
        ((patch_7, 0, +1), (patch_8, 0, -1), +1),
        ((patch_1, 1, +1), (patch_4, 1, -1), +1),
        ((patch_4, 1, +1), (patch_6, 1, -1), +1),
        ((patch_3, 1, +1), (patch_5, 1, -1), +1),
        ((patch_5, 1, +1), (patch_8, 1, -1), +1),
    )
    return Domain.join(patches, connectivity, name='square_8')


def build_square_9():
    """Build a square decomposed into a three-by-three patch grid."""
    logical_patch_1 = Square(
        'P1',
        bounds1=(0., np.pi / 3),
        bounds2=(0., np.pi / 3),
    )
    logical_patch_2 = Square(
        'P2',
        bounds1=(np.pi / 3, 2 * np.pi / 3),
        bounds2=(0., np.pi / 3),
    )
    logical_patch_3 = Square(
        'P3',
        bounds1=(2 * np.pi / 3, np.pi),
        bounds2=(0., np.pi / 3),
    )
    logical_patch_4 = Square(
        'P4',
        bounds1=(0., np.pi / 3),
        bounds2=(np.pi / 3, 2 * np.pi / 3),
    )
    logical_patch_5 = Square(
        'P5',
        bounds1=(2 * np.pi / 3, np.pi),
        bounds2=(np.pi / 3, 2 * np.pi / 3),
    )
    logical_patch_6 = Square(
        'P6',
        bounds1=(0., np.pi / 3),
        bounds2=(2 * np.pi / 3, np.pi),
    )
    logical_patch_7 = Square(
        'P7',
        bounds1=(np.pi / 3, 2 * np.pi / 3),
        bounds2=(2 * np.pi / 3, np.pi),
    )
    logical_patch_8 = Square(
        'P8',
        bounds1=(2 * np.pi / 3, np.pi),
        bounds2=(2 * np.pi / 3, np.pi),
    )
    logical_patch_9 = Square(
        'P9',
        bounds1=(np.pi / 3, 2 * np.pi / 3),
        bounds2=(np.pi / 3, 2 * np.pi / 3),
    )
    patch_1 = IdentityMapping('M1', dim=2)(logical_patch_1)
    patch_2 = IdentityMapping('M2', dim=2)(logical_patch_2)
    patch_3 = IdentityMapping('M3', dim=2)(logical_patch_3)
    patch_4 = IdentityMapping('M4', dim=2)(logical_patch_4)
    patch_5 = IdentityMapping('M5', dim=2)(logical_patch_5)
    patch_6 = IdentityMapping('M6', dim=2)(logical_patch_6)
    patch_7 = IdentityMapping('M7', dim=2)(logical_patch_7)
    patch_8 = IdentityMapping('M8', dim=2)(logical_patch_8)
    patch_9 = IdentityMapping('M9', dim=2)(logical_patch_9)

    patches = (
        patch_1, patch_2, patch_3, patch_4, patch_5,
        patch_6, patch_7, patch_8, patch_9,
    )
    connectivity = (
        ((patch_1, 0, +1), (patch_2, 0, -1), +1),
        ((patch_2, 0, +1), (patch_3, 0, -1), +1),
        ((patch_4, 0, +1), (patch_9, 0, -1), +1),
        ((patch_9, 0, +1), (patch_5, 0, -1), +1),
        ((patch_6, 0, +1), (patch_7, 0, -1), +1),
        ((patch_7, 0, +1), (patch_8, 0, -1), +1),
        ((patch_1, 1, +1), (patch_4, 1, -1), +1),
        ((patch_4, 1, +1), (patch_6, 1, -1), +1),
        ((patch_2, 1, +1), (patch_9, 1, -1), +1),
        ((patch_9, 1, +1), (patch_7, 1, -1), +1),
        ((patch_3, 1, +1), (patch_5, 1, -1), +1),
        ((patch_5, 1, +1), (patch_8, 1, -1), +1),
    )
    return Domain.join(patches, connectivity, name='square_9')


def build_annulus_3(r_min=None, r_max=None):
    """Build a complete annulus from three polar patches."""
    r_min = 0.5 if r_min is None else r_min
    r_max = 1.0 if r_max is None else r_max
    if not 0 < r_min < r_max:
        raise ValueError('annulus radii must satisfy 0 < r_min < r_max')

    logical_patch_1 = Square(
        'P1',
        bounds1=(r_min, r_max),
        bounds2=(0., np.pi / 2),
    )
    logical_patch_2 = Square(
        'P2',
        bounds1=(r_min, r_max),
        bounds2=(np.pi / 2, np.pi),
    )
    logical_patch_3 = Square(
        'P3',
        bounds1=(r_min, r_max),
        bounds2=(np.pi, 2 * np.pi),
    )
    mapping_1 = PolarMapping(
        'M1', dim=2, c1=0., c2=0., rmin=0., rmax=1.)
    mapping_2 = PolarMapping(
        'M2', dim=2, c1=0., c2=0., rmin=0., rmax=1.)
    mapping_3 = PolarMapping(
        'M3', dim=2, c1=0., c2=0., rmin=0., rmax=1.)
    patch_1 = mapping_1(logical_patch_1)
    patch_2 = mapping_2(logical_patch_2)
    patch_3 = mapping_3(logical_patch_3)

    patches = (patch_1, patch_2, patch_3)
    connectivity = (
        ((patch_1, 1, +1), (patch_2, 1, -1), +1),
        ((patch_2, 1, +1), (patch_3, 1, -1), +1),
        ((patch_3, 1, +1), (patch_1, 1, -1), +1),
    )
    return Domain.join(patches, connectivity, name='annulus_3')


def build_annulus_4(r_min=None, r_max=None):
    """Build a complete annulus from four quarter-annulus patches."""
    r_min = 0.5 if r_min is None else r_min
    r_max = 1.0 if r_max is None else r_max
    if not 0 < r_min < r_max:
        raise ValueError('annulus radii must satisfy 0 < r_min < r_max')

    logical_patch_1 = Square(
        'P1',
        bounds1=(r_min, r_max),
        bounds2=(0., np.pi / 2),
    )
    logical_patch_2 = Square(
        'P2',
        bounds1=(r_min, r_max),
        bounds2=(np.pi / 2, np.pi),
    )
    logical_patch_3 = Square(
        'P3',
        bounds1=(r_min, r_max),
        bounds2=(np.pi, 3 * np.pi / 2),
    )
    logical_patch_4 = Square(
        'P4',
        bounds1=(r_min, r_max),
        bounds2=(3 * np.pi / 2, 2 * np.pi),
    )
    mapping_1 = PolarMapping(
        'M1', dim=2, c1=0., c2=0., rmin=0., rmax=1.)
    mapping_2 = PolarMapping(
        'M2', dim=2, c1=0., c2=0., rmin=0., rmax=1.)
    mapping_3 = PolarMapping(
        'M3', dim=2, c1=0., c2=0., rmin=0., rmax=1.)
    mapping_4 = PolarMapping(
        'M4', dim=2, c1=0., c2=0., rmin=0., rmax=1.)
    patch_1 = mapping_1(logical_patch_1)
    patch_2 = mapping_2(logical_patch_2)
    patch_3 = mapping_3(logical_patch_3)
    patch_4 = mapping_4(logical_patch_4)

    patches = (patch_1, patch_2, patch_3, patch_4)
    connectivity = (
        ((patch_1, 1, +1), (patch_2, 1, -1), +1),
        ((patch_2, 1, +1), (patch_3, 1, -1), +1),
        ((patch_3, 1, +1), (patch_4, 1, -1), +1),
        ((patch_4, 1, +1), (patch_1, 1, -1), +1),
    )
    return Domain.join(patches, connectivity, name='annulus_4')


def build_curved_l_shape():
    """Build the three-patch curved L-shaped benchmark domain."""
    logical_patch_1 = Square(
        'P1',
        bounds1=(2, 3),
        bounds2=(0., np.pi / 8),
    )
    logical_patch_2 = Square(
        'P2',
        bounds1=(2, 3),
        bounds2=(np.pi / 8, np.pi / 4),
    )
    logical_patch_3 = Square(
        'P3',
        bounds1=(1, 2),
        bounds2=(np.pi / 8, np.pi / 4),
    )
    mapping_1 = PolarMapping(
        'M1', dim=2, c1=0., c2=0., rmin=0., rmax=1.)
    mapping_2 = PolarMapping(
        'M2', dim=2, c1=0., c2=0., rmin=0., rmax=1.)
    mapping_3 = PolarMapping(
        'M3', dim=2, c1=0., c2=0., rmin=0., rmax=1.)
    patch_1 = mapping_1(logical_patch_1)
    patch_2 = mapping_2(logical_patch_2)
    patch_3 = mapping_3(logical_patch_3)

    patches = (patch_1, patch_2, patch_3)
    connectivity = (
        ((patch_1, 1, +1), (patch_2, 1, -1), +1),
        ((patch_3, 0, +1), (patch_2, 0, -1), +1),
    )
    return Domain.join(patches, connectivity, name='curved_L_shape')


def build_pretzel(r_min=None, r_max=None):
    """Build the coarse eleven-patch pretzel domain."""
    r_min = 1 if r_min is None else r_min
    r_max = 2 if r_max is None else r_max
    if not 0 < r_min < r_max:
        raise ValueError('pretzel radii must satisfy 0 < r_min < r_max')

    h = r_max - r_min
    half_h = h / 2
    center_radius = h + (r_max + r_min) / 2

    logical_patch_1 = Square(
        'P1',
        bounds1=(r_min, r_max),
        bounds2=(0., np.pi / 2),
    )
    logical_patch_2 = Square(
        'P2',
        bounds1=(r_min, r_max),
        bounds2=(np.pi / 2, np.pi),
    )
    logical_patch_3 = Square(
        'P3',
        bounds1=(r_min, r_max),
        bounds2=(np.pi, 3 * np.pi / 2),
    )
    logical_patch_4 = Square(
        'P4',
        bounds1=(r_min, r_max),
        bounds2=(3 * np.pi / 2, 2 * np.pi),
    )
    mapping_1 = PolarMapping(
        'M1', dim=2, c1=h, c2=h, rmin=0., rmax=1.)
    mapping_2 = PolarMapping(
        'M2', dim=2, c1=-h, c2=h, rmin=0., rmax=1.)
    mapping_3 = PolarMapping(
        'M3', dim=2, c1=-h, c2=0, rmin=0., rmax=1.)
    mapping_4 = PolarMapping(
        'M4', dim=2, c1=h, c2=0, rmin=0., rmax=1.)
    patch_1 = mapping_1(logical_patch_1)
    patch_2 = mapping_2(logical_patch_2)
    patch_3 = mapping_3(logical_patch_3)
    patch_4 = mapping_4(logical_patch_4)

    logical_patch_5 = Square(
        'P5',
        bounds1=(-half_h, half_h),
        bounds2=(-h / 2, h / 2),
    )
    logical_patch_6 = Square(
        'P6',
        bounds1=(-half_h, half_h),
        bounds2=(-h / 2, h / 2),
    )
    logical_patch_7 = Square(
        'P7',
        bounds1=(-half_h, half_h),
        bounds2=(-h / 2, h / 2),
    )
    logical_patch_9 = Square(
        'P9',
        bounds1=(-half_h, half_h),
        bounds2=(-h, h),
    )
    logical_patch_12 = Square(
        'P12',
        bounds1=(-half_h, half_h),
        bounds2=(-h / 2, h / 2),
    )
    mapping_5 = AffineMapping(
        'M5', dim=2,
        c1=h / 2, c2=center_radius,
        a11=np.cos(np.pi / 2), a12=-np.sin(np.pi / 2),
        a21=np.sin(np.pi / 2), a22=np.cos(np.pi / 2),
    )
    mapping_6 = AffineMapping(
        'M6', dim=2,
        c1=-h / 2, c2=center_radius,
        a11=0, a12=1, a21=1, a22=0,
    )
    mapping_7 = AffineMapping(
        'M7', dim=2,
        c1=-center_radius, c2=h / 2,
        a11=np.cos(np.pi), a12=-np.sin(np.pi),
        a21=np.sin(np.pi), a22=np.cos(np.pi),
    )
    mapping_9 = AffineMapping(
        'M9', dim=2,
        c1=0, c2=h - center_radius,
        a11=np.cos(3 * np.pi / 2), a12=-np.sin(3 * np.pi / 2),
        a21=np.sin(3 * np.pi / 2), a22=np.cos(3 * np.pi / 2),
    )
    mapping_12 = AffineMapping(
        'M12', dim=2,
        c1=center_radius, c2=h / 2,
        a11=1, a12=0, a21=0, a22=-1,
    )
    patch_5 = mapping_5(logical_patch_5)
    patch_6 = mapping_6(logical_patch_6)
    patch_7 = mapping_7(logical_patch_7)
    patch_9 = mapping_9(logical_patch_9)
    patch_12 = mapping_12(logical_patch_12)

    logical_patch_13 = Square(
        'P13',
        bounds1=(3 * np.pi / 2, 2 * np.pi),
        bounds2=(r_min, r_max),
    )
    logical_patch_14 = Square(
        'P14',
        bounds1=(np.pi, 3 * np.pi / 2),
        bounds2=(r_min, r_max),
    )
    mapping_13 = TransposedPolarMapping(
        'M13', dim=2,
        c1=-r_min - h, c2=r_min + h, rmin=0., rmax=1.,
    )
    mapping_14 = TransposedPolarMapping(
        'M14', dim=2,
        c1=r_min + h, c2=r_min + h, rmin=0., rmax=1.,
    )
    patch_13 = mapping_13(logical_patch_13)
    patch_14 = mapping_14(logical_patch_14)

    patches = (
        patch_1, patch_2, patch_3, patch_4, patch_5, patch_6,
        patch_7, patch_9, patch_12, patch_13, patch_14,
    )
    connectivity = (
        ((patch_1, 1, +1), (patch_5, 1, -1), +1),
        ((patch_5, 1, +1), (patch_6, 1, +1), +1),
        ((patch_6, 1, -1), (patch_2, 1, -1), +1),
        ((patch_2, 1, +1), (patch_7, 1, -1), +1),
        ((patch_7, 1, +1), (patch_3, 1, -1), +1),
        ((patch_3, 1, +1), (patch_9, 1, -1), +1),
        ((patch_9, 1, +1), (patch_4, 1, -1), +1),
        ((patch_4, 1, +1), (patch_12, 1, +1), +1),
        ((patch_12, 1, -1), (patch_1, 1, -1), +1),
        ((patch_6, 0, -1), (patch_13, 0, +1), +1),
        ((patch_7, 0, -1), (patch_13, 0, -1), +1),
        ((patch_5, 0, -1), (patch_14, 0, -1), +1),
        ((patch_12, 0, -1), (patch_14, 0, +1), +1),
    )
    return Domain.join(patches, connectivity, name='pretzel')


def build_pretzel_annulus(r_min=None, r_max=None):
    """Build the nine-patch outer annular part of the pretzel domain."""
    r_min = 1 if r_min is None else r_min
    r_max = 2 if r_max is None else r_max
    if not 0 < r_min < r_max:
        raise ValueError('pretzel radii must satisfy 0 < r_min < r_max')

    h = r_max - r_min
    half_h = h / 2
    center_radius = h + (r_max + r_min) / 2

    logical_patch_1 = Square(
        'P1', bounds1=(r_min, r_max), bounds2=(0., np.pi / 2))
    logical_patch_2 = Square(
        'P2', bounds1=(r_min, r_max), bounds2=(np.pi / 2, np.pi))
    logical_patch_3 = Square(
        'P3', bounds1=(r_min, r_max),
        bounds2=(np.pi, 3 * np.pi / 2))
    logical_patch_4 = Square(
        'P4', bounds1=(r_min, r_max),
        bounds2=(3 * np.pi / 2, 2 * np.pi))
    mapping_1 = PolarMapping(
        'M1', dim=2, c1=h, c2=h, rmin=0., rmax=1.)
    mapping_2 = PolarMapping(
        'M2', dim=2, c1=-h, c2=h, rmin=0., rmax=1.)
    mapping_3 = PolarMapping(
        'M3', dim=2, c1=-h, c2=0, rmin=0., rmax=1.)
    mapping_4 = PolarMapping(
        'M4', dim=2, c1=h, c2=0, rmin=0., rmax=1.)
    patch_1 = mapping_1(logical_patch_1)
    patch_2 = mapping_2(logical_patch_2)
    patch_3 = mapping_3(logical_patch_3)
    patch_4 = mapping_4(logical_patch_4)

    logical_patch_5 = Square(
        'P5', bounds1=(-half_h, half_h), bounds2=(-h / 2, h / 2))
    logical_patch_6 = Square(
        'P6', bounds1=(-half_h, half_h), bounds2=(-h / 2, h / 2))
    logical_patch_7 = Square(
        'P7', bounds1=(-half_h, half_h), bounds2=(-h / 2, h / 2))
    logical_patch_9 = Square(
        'P9', bounds1=(-half_h, half_h), bounds2=(-h, h))
    logical_patch_12 = Square(
        'P12', bounds1=(-half_h, half_h), bounds2=(-h / 2, h / 2))
    mapping_5 = AffineMapping(
        'M5', dim=2,
        c1=h / 2, c2=center_radius,
        a11=np.cos(np.pi / 2), a12=-np.sin(np.pi / 2),
        a21=np.sin(np.pi / 2), a22=np.cos(np.pi / 2),
    )
    mapping_6 = AffineMapping(
        'M6', dim=2,
        c1=-h / 2, c2=center_radius,
        a11=0, a12=1, a21=1, a22=0,
    )
    mapping_7 = AffineMapping(
        'M7', dim=2,
        c1=-center_radius, c2=h / 2,
        a11=np.cos(np.pi), a12=-np.sin(np.pi),
        a21=np.sin(np.pi), a22=np.cos(np.pi),
    )
    mapping_9 = AffineMapping(
        'M9', dim=2,
        c1=0, c2=h - center_radius,
        a11=np.cos(3 * np.pi / 2), a12=-np.sin(3 * np.pi / 2),
        a21=np.sin(3 * np.pi / 2), a22=np.cos(3 * np.pi / 2),
    )
    mapping_12 = AffineMapping(
        'M12', dim=2,
        c1=center_radius, c2=h / 2,
        a11=1, a12=0, a21=0, a22=-1,
    )
    patch_5 = mapping_5(logical_patch_5)
    patch_6 = mapping_6(logical_patch_6)
    patch_7 = mapping_7(logical_patch_7)
    patch_9 = mapping_9(logical_patch_9)
    patch_12 = mapping_12(logical_patch_12)

    patches = (
        patch_1, patch_5, patch_6, patch_2, patch_7,
        patch_3, patch_9, patch_4, patch_12,
    )
    connectivity = (
        ((patch_1, 1, +1), (patch_5, 1, -1), +1),
        ((patch_5, 1, +1), (patch_6, 1, +1), +1),
        ((patch_6, 1, -1), (patch_2, 1, -1), +1),
        ((patch_2, 1, +1), (patch_7, 1, -1), +1),
        ((patch_7, 1, +1), (patch_3, 1, -1), +1),
        ((patch_3, 1, +1), (patch_9, 1, -1), +1),
        ((patch_9, 1, +1), (patch_4, 1, -1), +1),
        ((patch_4, 1, +1), (patch_12, 1, +1), +1),
        ((patch_12, 1, -1), (patch_1, 1, -1), +1),
    )
    return Domain.join(patches, connectivity, name='pretzel_annulus')


def build_pretzel_debug(r_min=None, r_max=None):
    """Build the two-patch pretzel geometry used for interface debugging."""
    r_min = 1 if r_min is None else r_min
    r_max = 2 if r_max is None else r_max
    if not 0 < r_min < r_max:
        raise ValueError('pretzel radii must satisfy 0 < r_min < r_max')

    h = r_max - r_min
    logical_patch_1 = Square(
        'P1',
        bounds1=(r_min, r_max),
        bounds2=(0., np.pi / 2),
    )
    logical_patch_10 = Square(
        'P10',
        bounds1=(r_min, r_max),
        bounds2=(np.pi / 2, np.pi),
    )
    mapping_1 = PolarMapping(
        'M1', dim=2, c1=h, c2=h, rmin=0., rmax=1.)
    mapping_10 = PolarMapping(
        'M10', dim=2, c1=h, c2=h, rmin=0., rmax=1.)
    patch_1 = mapping_1(logical_patch_1)
    patch_10 = mapping_10(logical_patch_10)

    patches = (patch_1, patch_10)
    connectivity = (
        ((patch_1, 1, +1), (patch_10, 1, -1), +1),
    )
    return Domain.join(patches, connectivity, name='pretzel_debug')


def build_pretzel_f(r_min=None, r_max=None):
    """Build the refined eighteen-patch pretzel domain."""
    r_min = 1 if r_min is None else r_min
    r_max = 2 if r_max is None else r_max
    if not 0 < r_min < r_max:
        raise ValueError('pretzel radii must satisfy 0 < r_min < r_max')

    h = r_max - r_min
    half_h = h / 2
    center_radius = h + (r_max + r_min) / 2

    logical_patch_1_1 = Square(
        'P1_1',
        bounds1=(r_min, r_max),
        bounds2=(0., np.pi / 4),
    )
    logical_patch_1_2 = Square(
        'P1_2',
        bounds1=(r_min, r_max),
        bounds2=(np.pi / 4, np.pi / 2),
    )
    logical_patch_2_1 = Square(
        'P2_1',
        bounds1=(r_min, r_max),
        bounds2=(np.pi / 2, 3 * np.pi / 4),
    )
    logical_patch_2_2 = Square(
        'P2_2',
        bounds1=(r_min, r_max),
        bounds2=(3 * np.pi / 4, np.pi),
    )
    logical_patch_3_1 = Square(
        'P3_1',
        bounds1=(r_min, r_max),
        bounds2=(np.pi, 5 * np.pi / 4),
    )
    logical_patch_3_2 = Square(
        'P3_2',
        bounds1=(r_min, r_max),
        bounds2=(5 * np.pi / 4, 3 * np.pi / 2),
    )
    logical_patch_4_1 = Square(
        'P4_1',
        bounds1=(r_min, r_max),
        bounds2=(3 * np.pi / 2, 7 * np.pi / 4),
    )
    logical_patch_4_2 = Square(
        'P4_2',
        bounds1=(r_min, r_max),
        bounds2=(7 * np.pi / 4, 2 * np.pi),
    )
    mapping_1_1 = PolarMapping(
        'M1_1', dim=2, c1=h, c2=h, rmin=0., rmax=1.)
    mapping_1_2 = PolarMapping(
        'M1_2', dim=2, c1=h, c2=h, rmin=0., rmax=1.)
    mapping_2_1 = PolarMapping(
        'M2_1', dim=2, c1=-h, c2=h, rmin=0., rmax=1.)
    mapping_2_2 = PolarMapping(
        'M2_2', dim=2, c1=-h, c2=h, rmin=0., rmax=1.)
    mapping_3_1 = PolarMapping(
        'M3_1', dim=2, c1=-h, c2=0, rmin=0., rmax=1.)
    mapping_3_2 = PolarMapping(
        'M3_2', dim=2, c1=-h, c2=0, rmin=0., rmax=1.)
    mapping_4_1 = PolarMapping(
        'M4_1', dim=2, c1=h, c2=0, rmin=0., rmax=1.)
    mapping_4_2 = PolarMapping(
        'M4_2', dim=2, c1=h, c2=0, rmin=0., rmax=1.)
    patch_1_1 = mapping_1_1(logical_patch_1_1)
    patch_1_2 = mapping_1_2(logical_patch_1_2)
    patch_2_1 = mapping_2_1(logical_patch_2_1)
    patch_2_2 = mapping_2_2(logical_patch_2_2)
    patch_3_1 = mapping_3_1(logical_patch_3_1)
    patch_3_2 = mapping_3_2(logical_patch_3_2)
    patch_4_1 = mapping_4_1(logical_patch_4_1)
    patch_4_2 = mapping_4_2(logical_patch_4_2)

    logical_patch_5 = Square(
        'P5',
        bounds1=(-half_h, half_h),
        bounds2=(-h / 2, h / 2),
    )
    logical_patch_6 = Square(
        'P6',
        bounds1=(-half_h, half_h),
        bounds2=(-h / 2, h / 2),
    )
    logical_patch_7 = Square(
        'P7',
        bounds1=(-half_h, half_h),
        bounds2=(-h / 2, h / 2),
    )
    logical_patch_9_1 = Square(
        'P9_1',
        bounds1=(-half_h, half_h),
        bounds2=(-h, 0),
    )
    logical_patch_9_2 = Square(
        'P9_2',
        bounds1=(-half_h, half_h),
        bounds2=(0, h),
    )
    logical_patch_12 = Square(
        'P12',
        bounds1=(-half_h, half_h),
        bounds2=(-h / 2, h / 2),
    )
    mapping_5 = AffineMapping(
        'M5', dim=2,
        c1=h / 2, c2=center_radius,
        a11=np.cos(np.pi / 2), a12=-np.sin(np.pi / 2),
        a21=np.sin(np.pi / 2), a22=np.cos(np.pi / 2),
    )
    mapping_6 = AffineMapping(
        'M6', dim=2,
        c1=-h / 2, c2=center_radius,
        a11=0, a12=1, a21=1, a22=0,
    )
    mapping_7 = AffineMapping(
        'M7', dim=2,
        c1=-center_radius, c2=h / 2,
        a11=np.cos(np.pi), a12=-np.sin(np.pi),
        a21=np.sin(np.pi), a22=np.cos(np.pi),
    )
    mapping_9_1 = AffineMapping(
        'M9_1', dim=2,
        c1=0, c2=h - center_radius,
        a11=np.cos(3 * np.pi / 2), a12=-np.sin(3 * np.pi / 2),
        a21=np.sin(3 * np.pi / 2), a22=np.cos(3 * np.pi / 2),
    )
    mapping_9_2 = AffineMapping(
        'M9_2', dim=2,
        c1=0, c2=h - center_radius,
        a11=np.cos(3 * np.pi / 2), a12=-np.sin(3 * np.pi / 2),
        a21=np.sin(3 * np.pi / 2), a22=np.cos(3 * np.pi / 2),
    )
    mapping_12 = AffineMapping(
        'M12', dim=2,
        c1=center_radius, c2=h / 2,
        a11=1, a12=0, a21=0, a22=-1,
    )
    patch_5 = mapping_5(logical_patch_5)
    patch_6 = mapping_6(logical_patch_6)
    patch_7 = mapping_7(logical_patch_7)
    patch_9_1 = mapping_9_1(logical_patch_9_1)
    patch_9_2 = mapping_9_2(logical_patch_9_2)
    patch_12 = mapping_12(logical_patch_12)

    logical_patch_13_1 = Square(
        'P13_1',
        bounds1=(3 * np.pi / 2, 7 * np.pi / 4),
        bounds2=(r_min, r_max),
    )
    logical_patch_13_2 = Square(
        'P13_2',
        bounds1=(7 * np.pi / 4, 2 * np.pi),
        bounds2=(r_min, r_max),
    )
    logical_patch_14_1 = Square(
        'P14_1',
        bounds1=(np.pi, 5 * np.pi / 4),
        bounds2=(r_min, r_max),
    )
    logical_patch_14_2 = Square(
        'P14_2',
        bounds1=(5 * np.pi / 4, 3 * np.pi / 2),
        bounds2=(r_min, r_max),
    )
    mapping_13_1 = TransposedPolarMapping(
        'M13_1', dim=2,
        c1=-r_min - h, c2=r_min + h, rmin=0., rmax=1.,
    )
    mapping_13_2 = TransposedPolarMapping(
        'M13_2', dim=2,
        c1=-r_min - h, c2=r_min + h, rmin=0., rmax=1.,
    )
    mapping_14_1 = TransposedPolarMapping(
        'M14_1', dim=2,
        c1=r_min + h, c2=r_min + h, rmin=0., rmax=1.,
    )
    mapping_14_2 = TransposedPolarMapping(
        'M14_2', dim=2,
        c1=r_min + h, c2=r_min + h, rmin=0., rmax=1.,
    )
    patch_13_1 = mapping_13_1(logical_patch_13_1)
    patch_13_2 = mapping_13_2(logical_patch_13_2)
    patch_14_1 = mapping_14_1(logical_patch_14_1)
    patch_14_2 = mapping_14_2(logical_patch_14_2)

    patches = (
        patch_1_1, patch_1_2, patch_2_1, patch_2_2,
        patch_3_1, patch_3_2, patch_4_1, patch_4_2,
        patch_5, patch_6, patch_7, patch_9_1, patch_9_2, patch_12,
        patch_13_1, patch_13_2, patch_14_1, patch_14_2,
    )
    connectivity = (
        ((patch_1_1, 1, +1), (patch_1_2, 1, -1), +1),
        ((patch_1_2, 1, +1), (patch_5, 1, -1), +1),
        ((patch_5, 1, +1), (patch_6, 1, +1), +1),
        ((patch_6, 1, -1), (patch_2_1, 1, -1), +1),
        ((patch_2_1, 1, +1), (patch_2_2, 1, -1), +1),
        ((patch_2_2, 1, +1), (patch_7, 1, -1), +1),
        ((patch_7, 1, +1), (patch_3_1, 1, -1), +1),
        ((patch_3_1, 1, +1), (patch_3_2, 1, -1), +1),
        ((patch_3_2, 1, +1), (patch_9_1, 1, -1), +1),
        ((patch_9_1, 1, +1), (patch_9_2, 1, -1), +1),
        ((patch_9_2, 1, +1), (patch_4_1, 1, -1), +1),
        ((patch_4_1, 1, +1), (patch_4_2, 1, -1), +1),
        ((patch_4_2, 1, +1), (patch_12, 1, +1), +1),
        ((patch_12, 1, -1), (patch_1_1, 1, -1), +1),
        ((patch_6, 0, -1), (patch_13_2, 0, +1), +1),
        ((patch_13_2, 0, -1), (patch_13_1, 0, +1), +1),
        ((patch_7, 0, -1), (patch_13_1, 0, -1), +1),
        ((patch_5, 0, -1), (patch_14_1, 0, -1), +1),
        ((patch_14_1, 0, +1), (patch_14_2, 0, -1), +1),
        ((patch_12, 0, -1), (patch_14_2, 0, +1), +1),
    )
    return Domain.join(patches, connectivity, name='pretzel_f')


def build_two_patch_twisted_strip():
    """Return two patches connected by interfaces of both 2D orientations."""
    logical_patch_a = Square('A', bounds1=(0, 1), bounds2=(0, 1))
    logical_patch_b = Square('B', bounds1=(0, 1), bounds2=(0, 1))

    patch_a = IdentityMapping('F_A', dim=2)(logical_patch_a)
    patch_b = AffineMapping(
        'F_B', dim=2,
        c1=1, c2=0,
        a11=1, a12=0,
        a21=0, a22=1,
    )(logical_patch_b)

    return Domain.join(
        patches=[patch_a, patch_b],
        interfaces=[
            # Touching sides preserve their tangential coordinates.
            ((patch_a, 0, +1), (patch_b, 0, -1), +1),
            # Outer sides are identified with reversed coordinates.
            ((patch_a, 0, -1), (patch_b, 0, +1), -1),
        ],
        name='two_patch_twisted_strip',
    )


def build_self_connected_patch():
    """Return a single square with a left-to-right self-interface."""
    logical_patch = Square('P', bounds1=(0, 1), bounds2=(0, 1))
    patch = IdentityMapping('F_P', dim=2)(logical_patch)

    return Domain.join(
        patches=[patch],
        interfaces=[((patch, 0, -1), (patch, 0, +1), +1)],
        name='self_connected_patch',
    )


def build_three_patch_shared_vertex():
    """Return three affine patches meeting at one interior vertex."""
    root_three_over_two = np.sqrt(3.0) / 2.0
    directions = (
        np.array((1.0, 0.0)),
        np.array((-0.5, root_three_over_two)),
        np.array((-0.5, -root_three_over_two)),
    )
    patches = []
    for patch_name, first, second in (
        ('A', directions[0], directions[1]),
        ('B', directions[1], directions[2]),
        ('C', directions[2], directions[0]),
    ):
        logical_patch = Square(
            patch_name, bounds1=(0, 1), bounds2=(0, 1))
        mapping = AffineMapping(
            f'F_{patch_name}', dim=2,
            c1=0, c2=0,
            a11=first[0], a12=second[0],
            a21=first[1], a22=second[1],
        )
        patches.append(mapping(logical_patch))

    patch_a, patch_b, patch_c = patches
    return Domain.join(
        patches=patches,
        interfaces=[
            ((patch_a, 0, -1), (patch_b, 1, -1), +1),
            ((patch_b, 0, -1), (patch_c, 1, -1), +1),
            ((patch_c, 0, -1), (patch_a, 1, -1), +1),
        ],
        name='three_patch_shared_vertex',
    )


_DOMAIN_BUILDERS_2D = {
    'square_2': build_square_2,
    'square_4': build_square_4,
    'square_6': build_square_6,
    'square_8': build_square_8,
    'square_9': build_square_9,
    'annulus_3': build_annulus_3,
    'annulus_4': build_annulus_4,
    'curved_L_shape': build_curved_l_shape,
    'pretzel': build_pretzel,
    'pretzel_f': build_pretzel_f,
    'pretzel_annulus': build_pretzel_annulus,
    'pretzel_debug': build_pretzel_debug,
    'two_patch_twisted_strip': build_two_patch_twisted_strip,
    'self_connected_patch': build_self_connected_patch,
    'three_patch_shared_vertex': build_three_patch_shared_vertex,
}

_RADIAL_DOMAINS_2D = {
    'annulus_3',
    'annulus_4',
    'pretzel',
    'pretzel_f',
    'pretzel_annulus',
    'pretzel_debug',
}


def build_multipatch_domain_2d(domain_name='square_2', r_min=None, r_max=None):
    """Build a named 2D multipatch domain from the gallery.

    Parameters
    ----------
    domain_name : str
        Name registered in the 2D multipatch gallery.
    r_min, r_max : float, optional
        Inner and outer radii for annulus and pretzel domains.
    """
    try:
        builder = _DOMAIN_BUILDERS_2D[domain_name]
    except (KeyError, TypeError) as error:
        choices = ', '.join(_DOMAIN_BUILDERS_2D)
        raise ValueError(
            f'unknown 2D multipatch domain {domain_name!r}; '
            f'choose from {choices}') from error

    if domain_name in _RADIAL_DOMAINS_2D:
        return builder(r_min=r_min, r_max=r_max)
    return builder()


def build_cartesian_multipatch_domain_2d(
        ncells, log_interval_x, log_interval_y, mapping='identity'):
    """Create a 2D multipatch domain from a rectangular patch layout.

    Parameters
    ----------
    ncells : array-like
        Two-dimensional patch layout. Non-``None`` entries create patches;
        their numerical values are ignored by this symbolic builder and may
        subsequently be used as per-patch cell counts. Row zero is the top
        row of the domain.
    log_interval_x : tuple
        The interval in the x direction in the logical domain.
    log_interval_y : tuple
        The interval in the y direction in the logical domain.
    mapping : str
        The type of mapping to use. Can be ``identity`` or ``polar``.

    Returns
    -------
    domain : Domain
        The symbolic multipatch domain.
    """
    layout = np.asarray(ncells, dtype=object)
    if layout.ndim != 2:
        raise ValueError('ncells must be a two-dimensional patch layout')

    nrows, ncols = layout.shape
    if nrows == 0 or ncols == 0:
        raise ValueError('ncells must contain at least one layout position')

    try:
        ax, bx = log_interval_x
        ay, by = log_interval_y
    except (TypeError, ValueError) as error:
        raise ValueError(
            'logical intervals must each contain exactly two bounds') from error

    mapping_types = {
        'identity': IdentityMapping,
        'polar': PolarMapping,
    }
    try:
        mapping_type = mapping_types[mapping]
    except (KeyError, TypeError) as error:
        raise ValueError(
            "mapping must be either 'identity' or 'polar'") from error

    patches_by_position = {}
    patches = []
    for row in range(nrows):
        bounds2 = (
            by - (row + 1) / nrows * (by - ay),
            by - row / nrows * (by - ay),
        )
        for column in range(ncols):
            if layout[row, column] is None:
                continue

            bounds1 = (
                ax + column / ncols * (bx - ax),
                ax + (column + 1) / ncols * (bx - ax),
            )
            logical_patch = Square(
                f'Log_{row}_{column}',
                bounds1=bounds1,
                bounds2=bounds2,
            )
            mapping_kwargs = {}
            if mapping == 'polar':
                mapping_kwargs = dict(c1=0., c2=0., rmin=0., rmax=1.)
            patch_mapping = mapping_type(
                f'M_{row}_{column}', dim=2, **mapping_kwargs)
            patch = patch_mapping(logical_patch)
            patches_by_position[row, column] = patch
            patches.append(patch)

    if not patches:
        raise ValueError('ncells must contain at least one non-None entry')

    connectivity = []
    for row in range(nrows):
        for column in range(ncols - 1):
            left = patches_by_position.get((row, column))
            right = patches_by_position.get((row, column + 1))
            if left is not None and right is not None:
                connectivity.append(
                    ((left, 0, +1), (right, 0, -1), +1))

    for row in range(nrows - 1):
        for column in range(ncols):
            top = patches_by_position.get((row, column))
            bottom = patches_by_position.get((row + 1, column))
            if top is not None and bottom is not None:
                connectivity.append(
                    ((top, 1, -1), (bottom, 1, +1), +1))

    return Domain.join(patches, connectivity, name='domain')


# =============================================================================
# 3D domains
# =============================================================================

def build_two_cube_3d():
    """Build two mapped cubes joined across differently parameterized faces."""
    logical_patch_a = Cube('A')
    logical_patch_b = Cube('B')

    mapping_a = AffineMapping(
        'F_A', dim=3,
        c1=0.0, c2=0.0, c3=0.0,
        a11=1.0, a12=0.25, a13=0.0,
        a21=0.2, a22=1.00, a23=0.15,
        a31=0.0, a32=0.10, a33=1.0,
    )
    mapping_b = AffineMapping(
        'F_B', dim=3,
        c1=1.0, c2=0.35, c3=1.0,
        a11=0.00, a12=1.0, a13=0.25,
        a21=-0.15, a22=0.2, a23=1.00,
        a31=-1.00, a32=0.0, a33=0.10,
    )
    patch_a = mapping_a(logical_patch_a)
    patch_b = mapping_b(logical_patch_b)

    patches = (patch_a, patch_b)
    connectivity = (
        (
            (patch_a, 0, +1),
            (patch_b, 1, -1),
            (-1, +1, -1),
        ),
    )
    return Domain.join(
        patches,
        connectivity,
        name='mapped_two_cube_3d',
    )


def build_torus_3d(
        major_radius=2.0, minor_bounds=(0.35, 0.80), *, hollow=True,
        toroidal_angle=2.0 * np.pi, close_torus=None):
    """Build a four-patch torus with two angular cuts per direction.

    Parameters
    ----------
    major_radius : float, default=2.0
        Distance from the centre of the tube to the axis of revolution.
    minor_bounds : tuple[float, float], default=(0.35, 0.80)
        Inner and outer minor radii used by the regular hollow variant. The
        lower value must be positive and the upper value must be smaller than
        ``major_radius``.
    hollow : bool, default=True
        Keep the positive inner minor radius. If false, replace it by zero to
        fill the tube. The resulting mapping is singular on the centreline.
    toroidal_angle : float, default=2*pi
        Total angle swept around the axis of revolution, in radians. It must
        be greater than zero and at most ``2*pi``. The two toroidal patch
        sectors each cover half of this angle.
    close_torus : bool or None, default=None
        Join the two toroidal end faces. ``None`` closes a full ``2*pi`` sweep
        and leaves every smaller sweep open. Explicit closure is valid only
        for a full sweep; ``False`` may also leave a full sweep cut open.
    """
    inner_radius, outer_radius = map(float, minor_bounds)
    major_radius = float(major_radius)
    if not 0.0 < inner_radius < outer_radius < major_radius:
        raise ValueError(
            'expected 0 < inner minor radius < outer minor radius '
            '< major radius')
    if not isinstance(hollow, bool):
        raise TypeError('hollow must be a bool')
    if not hollow:
        inner_radius = 0.0

    full_angle = 2.0 * np.pi
    toroidal_angle = float(toroidal_angle)
    if not 0.0 < toroidal_angle <= full_angle or not np.isfinite(
            toroidal_angle):
        raise ValueError('toroidal_angle must satisfy 0 < angle <= 2*pi')
    is_full_sweep = bool(np.isclose(toroidal_angle, full_angle))
    if is_full_sweep:
        toroidal_angle = full_angle
    if close_torus is None:
        close_torus = is_full_sweep
    elif not isinstance(close_torus, bool):
        raise TypeError('close_torus must be a bool or None')
    if close_torus and not is_full_sweep:
        raise ValueError('only a full 2*pi sweep can be closed')

    theta_bounds = (0.0, np.pi, full_angle)
    phi_bounds = (0.0, 0.5 * toroidal_angle, toroidal_angle)
    patch_grid = {}
    patches = []

    for theta_sector in range(2):
        for phi_sector in range(2):
            patch_name = f'P{theta_sector}{phi_sector}'
            logical_patch = Cube(
                patch_name,
                bounds1=(inner_radius, outer_radius),
                bounds2=theta_bounds[theta_sector:theta_sector + 2],
                bounds3=phi_bounds[phi_sector:phi_sector + 2],
            )
            mapping = TorusMapping(
                f'F_{patch_name}', dim=3, R0=major_radius)
            patch = mapping(logical_patch)
            patch_grid[theta_sector, phi_sector] = patch
            patches.append(patch)

    identity = (+1, +1, +1)
    connectivity = []
    for phi_sector in range(2):
        lower = patch_grid[0, phi_sector]
        upper = patch_grid[1, phi_sector]
        connectivity.extend((
            ((lower, 1, +1), (upper, 1, -1), identity),
            ((upper, 1, +1), (lower, 1, -1), identity),
        ))

    for theta_sector in range(2):
        lower = patch_grid[theta_sector, 0]
        upper = patch_grid[theta_sector, 1]
        connectivity.append(
            ((lower, 2, +1), (upper, 2, -1), identity))
        if close_torus:
            connectivity.append(
                ((upper, 2, +1), (lower, 2, -1), identity))

    shape = 'hollow' if hollow else 'solid'
    closure = '' if close_torus else 'open_'
    return Domain.join(
        patches,
        connectivity,
        name=f'{closure}{shape}_torus_3d',
    )


def build_oriented_two_cube():
    """Return two cubes with a swapped and reversed face parameterization."""
    logical_patch_a = Cube('A')
    logical_patch_b = Cube('B')

    mapping_a = IdentityMapping('F_A', dim=3)
    mapping_b = AffineMapping(
        'F_B', dim=3,
        c1=1, c2=0, c3=1,
        a11=0, a12=1, a13=0,
        a21=0, a22=0, a23=1,
        a31=-1, a32=0, a33=0,
    )
    patch_a = mapping_a(logical_patch_a)
    patch_b = mapping_b(logical_patch_b)

    return Domain.join(
        patches=[patch_a, patch_b],
        interfaces=[(
            (patch_a, 0, +1),
            (patch_b, 1, -1),
            (-1, +1, -1),
        )],
        name='oriented_two_cube',
    )


def build_three_patch_shared_edge():
    """Return three affine hexahedra meeting at one interior edge."""
    root_three_over_two = np.sqrt(3.0) / 2.0
    directions = (
        np.array((1.0, 0.0)),
        np.array((-0.5, root_three_over_two)),
        np.array((-0.5, -root_three_over_two)),
    )

    patches = []
    for patch_name, first, second in (
        ('H0', directions[0], directions[1]),
        ('H1', directions[1], directions[2]),
        ('H2', directions[2], directions[0]),
    ):
        logical_patch = Cube(patch_name)
        mapping = AffineMapping(
            f'G_{patch_name}', dim=3,
            c1=0, c2=0, c3=0,
            a11=first[0], a12=second[0], a13=0,
            a21=first[1], a22=second[1], a23=0,
            a31=0, a32=0, a33=1,
        )
        patches.append(mapping(logical_patch))

    patch_0, patch_1, patch_2 = patches
    return Domain.join(
        patches=patches,
        interfaces=[
            ((patch_0, 0, -1), (patch_1, 1, -1), (+1, +1, +1)),
            ((patch_1, 0, -1), (patch_2, 1, -1), (+1, +1, +1)),
            ((patch_2, 0, -1), (patch_0, 1, -1), (+1, +1, +1)),
        ],
        name='three_patch_shared_edge',
    )


_DOMAIN_BUILDERS_3D = {
    'two_cube': build_two_cube_3d,
    'torus': build_torus_3d,
    'oriented_two_cube': build_oriented_two_cube,
    'three_patch_shared_edge': build_three_patch_shared_edge,
}


def build_multipatch_domain_3d(domain_name='two_cube', **kwargs):
    """Build a named 3D multipatch domain from the gallery.

    Parameters
    ----------
    domain_name : str
        Name registered in the 3D multipatch gallery.
    kwargs : object
        Optional arguments passed to the selected domain builder.
    """
    try:
        builder = _DOMAIN_BUILDERS_3D[domain_name]
    except (KeyError, TypeError) as error:
        choices = ', '.join(_DOMAIN_BUILDERS_3D)
        raise ValueError(
            f'unknown 3D multipatch domain {domain_name!r}; '
            f'choose from {choices}') from error

    return builder(**kwargs)


# =============================================================================
# Plotting
# =============================================================================

def plot_multipatch_domain(
        domain_name='torus', *, output=None, show=True, topology=False,
        builder_options=None):
    """Build, print, and plot a named 2D or 3D gallery domain.

    Parameters
    ----------
    domain_name : str
        Name registered in the 2D or 3D multipatch gallery.
    output : path-like or None, default=None
        Optional path of the generated PNG or PDF figure.
    show : bool, default=True
        Open an interactive Matplotlib window.
    topology : bool, default=False
        Add patch, interface, orientation, and vertex annotations.
    builder_options : dict or None, default=None
        Optional arguments passed to the selected domain builder.

    Returns
    -------
    domain : Domain
        The selected symbolic multipatch domain.
    figure : matplotlib.figure.Figure
        The generated Matplotlib figure.
    """
    options = {} if builder_options is None else dict(builder_options)
    if domain_name in _DOMAIN_BUILDERS_2D:
        domain = build_multipatch_domain_2d(domain_name, **options)
    elif domain_name in _DOMAIN_BUILDERS_3D:
        domain = build_multipatch_domain_3d(domain_name, **options)
    else:
        choices = ', '.join((*_DOMAIN_BUILDERS_2D, *_DOMAIN_BUILDERS_3D))
        raise ValueError(
            f'unknown multipatch domain {domain_name!r}; '
            f'choose from {choices}')

    import matplotlib.pyplot as plt
    from sympde.utilities import plot_domain, print_topology

    print_topology(domain)
    figure = plot_domain(
        domain,
        draw=False,
        isolines=not topology,
        topology=topology,
        interface_labels=topology,
        vertex_labels=topology,
    )
    if not topology:
        figure.axes[0].set_title(str(domain.name).replace('_', ' '))

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=180, bbox_inches='tight')
        print(f'Saved plot to {output}')

    if show:
        plt.show()
    else:
        plt.close(figure)

    return domain, figure


def main(argv=None):
    """Plot a gallery domain selected from the command line."""
    domain_names = (*_DOMAIN_BUILDERS_2D, *_DOMAIN_BUILDERS_3D)
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        'domain_name', nargs='?', default='torus', choices=domain_names,
        help='gallery domain to build and plot (default: torus)')
    parser.add_argument(
        '--output', type=Path,
        help='optional path of the generated PNG or PDF figure')
    parser.add_argument(
        '--no-show', action='store_true',
        help='create the plot without opening an interactive window')
    parser.add_argument(
        '--topology', action='store_true',
        help='add patch, interface, orientation, and vertex annotations')
    arguments = parser.parse_args(argv)

    return plot_multipatch_domain(
        arguments.domain_name,
        output=arguments.output,
        show=not arguments.no_show,
        topology=arguments.topology,
    )


if __name__ == '__main__':
    main()
