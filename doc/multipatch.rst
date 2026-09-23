Multipatch domains
******************

This chapter describes the symbolic multipatch topology in SymPDE.  No mesh,
spline space, degrees of freedom, or MPI decomposition is created here.  A
SymPDE multipatch domain answers topological questions such as:

* which tensor-product patches constitute the domain;
* which boundary faces are identified;
* how the tangential coordinates of the two faces correspond;
* which faces remain on the exterior boundary; and
* which patch-local corners represent the same topological vertex.

Psydac consumes this information later when it creates a discrete geometry.
Keeping the topology independent of a particular discretization is important:
the same domain may be discretized with different cell counts, degrees, knot
vectors, and process counts.

The stable collection views described below are additive, and ordinary
single-patch construction with ``Line``, ``Square``, ``Cube``, mappings,
``domain.interior``, and ``domain.boundary`` is unchanged.  The multipatch
interface API and serialized interface format intentionally replace older,
ambiguous representations; see `Migration from the previous multipatch API`_.

The object model
================

A multipatch domain has four important levels::

   logical NCube patches              Square('A'), Cube('B'), ...
           |
           | optional patch mappings
           v
   physical patches                   F_A(A), F_B(B), ...
           |
           | Domain.join(...)
           v
   Domain
      |-- patches                     stable tuple of patch interiors
      |-- boundary                    faces not used by an interface
      |-- interface_map               name -> Interface
      |-- exterior_sides              stable tuple of exterior faces
      `-- shared_vertices             stable tuple of vertex classes

These descriptive properties are views of the existing symbolic objects.
``domain.interior``, ``domain.subdomains``, ``domain.connectivity``,
``domain.boundary``, and ``domain.corners`` retain their original behavior.

The logical dimension is denoted by ``dim`` below.  A patch face has dimension
``dim - 1``.  Logical axes are zero-based: a square has axes 0 and 1, and a
cube has axes 0, 1, and 2.

Creating patches
================

The usual logical patches are ``Line``, ``Square``, and ``Cube``:

.. code-block:: python

   from sympde.topology import Line, Square, Cube

   interval = Line('I', bounds=(-1, 2))
   square = Square('A', bounds1=(0, 2), bounds2=(-1, 1))
   cube = Cube(
       'B',
       bounds1=(0, 1),
       bounds2=(0, 2),
       bounds3=(-1, 1),
   )

The parameters have the following meanings.

``name``
   A unique patch name.  Names are used for inspection, serialization, and
   patch-specific discretization parameters.

``bounds``, ``bounds1``, ``bounds2``, ``bounds3``
   The lower and upper logical coordinate of each tensor-product direction.
   Every lower bound must be strictly smaller than its upper bound.  Bounds
   describe the parameter domain; they do not place a mapped patch in physical
   space.

``NCube(name, dim, min_coords, max_coords)`` is the dimension-independent
constructor.  ``Line``, ``Square``, and ``Cube`` are its convenient 1D, 2D,
and 3D forms.  Multipatch geometry is developed and tested primarily for these
three dimensions.

Logical and physical patches
----------------------------

A mapping is applied to a logical patch by calling it:

.. code-block:: python

   from sympde.topology import AffineMapping, IdentityMapping, Square

   logical_a = Square('A')
   logical_b = Square('B')

   patch_a = IdentityMapping('F_A', dim=2)(logical_a)
   patch_b = AffineMapping(
       'F_B', dim=2,
       c1=1, c2=0,
       a11=1, a12=0,
       a21=0, a22=1,
   )(logical_b)

For ``Mapping(name, dim=dim)``, ``dim`` sets both the logical and physical
dimension.  More generally, ``ldim`` and ``pdim`` may be supplied by mapping
classes which support embedded geometries.  Parameters such as ``c1`` and
``a11`` are specific to the mapping class.  An ``AffineMapping`` in 2D is

.. math::

   F(\eta_1,\eta_2) =
   \begin{pmatrix}c_1\\c_2\end{pmatrix} +
   \begin{pmatrix}a_{11}&a_{12}\\a_{21}&a_{22}\end{pmatrix}
   \begin{pmatrix}\eta_1\\\eta_2\end{pmatrix}.

If every patch passed to ``Domain.join`` is mapped, the joined domain receives
a ``MultiPatchMapping`` and a corresponding ``logical_domain``.  Mixing mapped
and unmapped patches is representable, but no aggregate multipatch mapping is
created.  Downstream geometry and variational code generally expects the
patches to be consistently all mapped or all logical.

A mapping may also be applied to an already joined logical domain.  In that
case its interfaces are mapped together with its interiors and boundaries.
Each mapped interface retains the original interface as its ``logical_domain``
and receives an ``InterfaceMapping`` for interface integrals.

Domain gallery
==============

Reusable symbolic examples live in ``sympde.topology.multipatch_gallery``.
They have no dependency on Psydac or on a discretization:

.. code-block:: python

   from sympde.topology.multipatch_gallery import (
       build_cartesian_multipatch_domain_2d,
       build_multipatch_domain_2d,
   )

   annulus = build_multipatch_domain_2d('annulus_4')
   square_with_a_hole = build_multipatch_domain_2d('square_8')

``build_multipatch_domain_2d`` also provides ``square_2``, ``square_4``,
``square_6``, ``square_9``, ``annulus_3``, ``curved_L_shape``, ``pretzel``,
``pretzel_f``, ``pretzel_annulus``, and ``pretzel_debug``.
Each named geometry has a dedicated public builder in the same module; the
name-based function is only a convenience dispatcher.

``build_cartesian_multipatch_domain_2d`` creates a rectangular grid from a 2D
layout. Every non-``None`` entry creates one patch and a ``None`` entry leaves
a hole. The values themselves are not interpreted by SymPDE, so the same array
may later carry patch-specific cell counts for a discretizer. Both identity
and polar patch mappings are available.

The gallery also contains two 3D domains. ``build_two_patch_3d`` constructs
two sheared cubes with a cross-axis interface, and ``build_torus_2x2_3d``
constructs the four-patch hollow or solid toroidal volume.
``build_multipatch_domain_3d`` selects these builders by the names
``two_patch`` and ``torus_2x2`` and forwards optional arguments to the
selected builder.

The gallery module can also plot any registered 2D or 3D domain directly.
For example, the following commands open an ordinary annulus plot and save an
annotated two-cube topology plot, respectively:

.. code-block:: console

   python -m sympde.topology.multipatch_gallery annulus_4
   python sympde/topology/multipatch_gallery.py two_patch \
       --topology --no-show --output two_patch.png

The corresponding Python entry point is ``plot_multipatch_domain``. Its
``builder_options`` argument forwards domain-specific parameters to the
selected 2D or 3D builder.

Joining patches
===============

The public constructor is

.. code-block:: python

   Domain.join(patches, interfaces=interfaces, name=name)

``patches``
   A non-empty ordered list or tuple of ``Domain`` objects.  All patches must
   have the same logical dimension.  The order defines the integer patch
   indices accepted in the interface specifications.  Patch objects, physical
   interior names, and logical names used by serialization must all be unique.

``interfaces``
   A list or tuple of interface descriptions.  Each description is the tuple
   ``(minus, plus, orientation)`` or
   ``(minus, plus, orientation, name)``.  The optional explicit name is
   preserved exactly and must be unique.  The whole collection is a sequence
   because several distinct interfaces may connect the same pair of patches.
   ``connectivity`` is an equivalent positional/keyword spelling.
   Omitting both arguments is equivalent to providing an empty tuple.

``name``
   The name of the resulting domain.

If one patch and an empty connectivity are supplied, that patch is returned
unchanged.  A one-patch domain with non-empty connectivity follows the normal
join path and can therefore contain self-interfaces.

Describing an interface side
----------------------------

Each side is the triple ``(patch, axis, ext)``:

``patch``
   Either the patch object itself or its integer position in ``patches``.

``axis``
   The logical direction normal to the face, in ``range(dim)``.

``ext``
   ``-1`` selects the minimum-coordinate face and ``+1`` the
   maximum-coordinate face.

For example:

.. code-block:: python

   left = (A, 0, -1)
   top = (A, 1, +1)

An interface combines two side triples with their orientation:

.. code-block:: python

   interface = ((A, 0, +1), (B, 1, -1), -1)

``Domain.join`` resolves each side to an existing symbolic ``Boundary``
object.

``minus`` and ``plus`` are ordered labels.  They do not mean that one patch is
geometrically smaller, coarser, or located on a particular side in physical
space.  The ordering establishes the direction in which the orientation is
read and remains the ordering stored by the resulting ``Interface``.

Interface orientation
=====================

The public orientation value depends on the patch dimension:

``1D``
   ``None``.  An endpoint has no tangential coordinate.

``2D``
   One sign, ``+1`` or ``-1``.  It says whether the coordinates along the two
   edges increase in the same or opposite directions.

``3D``
   ``(flag, sign1, sign2)``, where all three entries are ``+1`` or ``-1``.
   ``flag=+1`` preserves the order of the two tangential coordinates;
   ``flag=-1`` exchanges their order.  ``sign1`` and ``sign2`` then give the
   directions of the first and second resulting correspondences.

Thus the common orientations are:

.. code-block:: python

   identity_2d = +1
   reversed_2d = -1

   identity_3d = (+1, +1, +1)
   swap_3d = (-1, +1, +1)
   swap_and_flip_3d = (-1, +1, -1)

The ordered tangential axes are obtained by removing the face-normal axis from
``range(dim)`` without otherwise changing the order.  The most convenient
way to inspect the decoded result is ``interface.axis_map``.  It returns
``(minus_axis, plus_axis, direction)`` for every tangential direction.

For example, consider a 3D interface whose minus normal is axis 0 and whose
plus normal is axis 1:

.. code-block:: python

   from sympde.topology import Cube, Domain

   patch_a = Cube('A')
   patch_b = Cube('B')
   domain = Domain.join(
       patches=[patch_a, patch_b],
       interfaces=[(
           (0, 0, +1),
           (1, 1, -1),
           (-1, +1, -1),
       )],
       name='Omega',
   )
   interface, = domain.interface_map.values()
   assert interface.orientation == (-1, +1, -1)
   assert interface.axis_map == ((1, 2, 1), (2, 0, -1))

The first tuple says minus axis 1 maps to plus axis 2 in the same direction.
The second says minus axis 2 maps to plus axis 0 in the opposite direction.
Normal axes are intentionally absent from the orientation value: they already
belong to the two side triples.  The orientation is decoded against each
face's own normal axis.  Therefore the two normal axes need not be equal.

An ``Interface`` deliberately has no single ``axis`` property because its two
faces may have different normal axes.  Read ``interface.minus.axis`` and
``interface.plus.axis`` explicitly.  New multipatch code may use the clearer
aliases ``interface.minus_side.normal_axis``,
``interface.plus_side.normal_axis``, and ``interface.axis_map``.

Identity orientation means that tangential positions retain their order and
all signs are positive.  It does not mean that the two normal axes must be
equal, their extremities must differ, or the physical mappings must be
identical.

A complete 2D example
=====================

This example joins two mapped squares twice.  The touching seam preserves its
tangential coordinate, while the outer seam reverses it:

.. code-block:: python

   from sympde.topology import (
       AffineMapping, Domain, IdentityMapping, Square,
   )

   A = IdentityMapping('F_A', dim=2)(Square('A'))
   B = AffineMapping(
       'F_B', dim=2,
       c1=1, c2=0,
       a11=1, a12=0,
       a21=0, a22=1,
   )(Square('B'))

   omega = Domain.join(
       patches=[A, B],
       interfaces=[
           ((A, 0, +1), (B, 0, -1), +1),
           ((A, 0, -1), (B, 0, +1), -1),
       ],
       name='twisted_strip',
   )

This example illustrates two independent facts:

* the same ordered patch pair can have more than one interface; and
* the orientation belongs to each interface, not to the patch pair.

Interface names are generated from the side-domain names.  If that base name
is repeated, later occurrences receive ``#2``, ``#3``, and so on.  The names
are checked against all previously generated names, including names which
already contain ``#`` or ``|``, so one interface cannot overwrite another.
They are stable keys in ``omega.interface_map``.  ``omega.connectivity`` is
the unchanged symbolic storage object behind that descriptive view.

Self-interfaces
===============

Two distinct faces of one patch may be identified:

.. code-block:: python

   P = Square('P')
   periodic_x = Domain.join(
       patches=[P],
       interfaces=[((P, 0, -1), (P, 0, +1), +1)],
       name='periodic_x',
   )

Both ``interface.minus.domain`` and ``interface.plus.domain`` are ``P``.  The
two boundary objects remain distinct because their extremities differ.  A
self-interface is part of the topology; whether a particular discretization
consumer realizes it as periodic basis functions or as an explicit constraint
is a separate question.

Shared vertices and corners
===========================

Joining faces also identifies their corners.  SymPDE computes this transitively:
if an interface identifies a corner of A with a corner of B and another
identifies that B corner with a corner of C, all three belong to one
shared vertex.  ``SharedVertex`` is the descriptive alias of the historical
``CornerInterface`` class, while ``PatchVertex`` aliases the patch-local
``CornerBoundary`` class.  The historical names remain canonical for SymPy
printing and serialization.

For each interface, SymPDE enumerates the ``2**(dim - 1)`` corners of its minus
face, maps their tangential extremities with
``interface.axis_map``, and unions the resulting minus and plus
``PatchVertex`` objects.  This supports, for
example, three or more 2D patches meeting at one vertex and 3D faces with
permuted axes.

Use ``domain.interface_map`` and ``domain.shared_vertices`` for stable
multipatch collection views:

.. code-block:: python

   for name, interface in omega.interface_map.items():
       print(name)
       print(interface.minus_side.patch,
             interface.minus_side.normal_axis,
             interface.minus_side.ext)
       print(interface.plus_side.patch,
             interface.plus_side.normal_axis,
             interface.plus_side.ext)
       print(interface.orientation)
       print(interface.axis_map)

   for vertex in omega.shared_vertices:
       for local_vertex in vertex.corners:
           print(local_vertex.domain, local_vertex.coordinates)

``PatchVertex.coordinates`` is a tuple of zeros and ones, one per logical
axis.  It identifies the minimum or maximum vertex in patch-local coordinates.
``domain.shared_vertices`` contains shared vertices induced by connectivity;
it is not a list of every unshared exterior vertex.  The symbolic
``domain.corners`` property remains available with its original
``None``/object/``Union`` cardinality behavior.

Inspecting and plotting 2D and 3D topology
==========================================

The plotting utilities can turn the symbolic connectivity into structured
vertex data, a textual summary, or an annotated figure.  No mesh or Psydac
discretization is created:

.. code-block:: python

   from sympde.utilities import (
       collect_topology_vertices,
       plot_domain,
       print_topology,
   )

   vertices = collect_topology_vertices(omega)
   for vertex in vertices:
       print(vertex.index, vertex.patches, vertex.is_boundary)
       for incidence in vertex.incidences:
           print(
               incidence.patch_index,
               incidence.patch_name,
               incidence.logical_corner,
           )

   print_topology(omega)

   figure = plot_domain(
       omega,
       draw=False,
       refinement=60,
       isolines=True,
       topology=True,
       vertex_labels=True,
   )
   figure.savefig('omega-topology.png', dpi=180, bbox_inches='tight')

``collect_topology_vertices`` differs deliberately from ``domain.corners``:
it begins with every corner of every patch and merges the shared-corner
equivalence classes.  Its result therefore includes unshared exterior
vertices as well as interior or multiply represented vertices.  Patch names
and integer patch indices are both recorded.

In a 2D annotated topology plot, the minus interface edge is solid, the plus
edge is dashed, and arrows point along each side's increasing native
tangential coordinate. Interface colors pair the two sides, while the text
legend reports the signed ``axis_map``. If identified sides are separated in
the drawing, a dotted curve connects their midpoints; this makes cut-open
periodic or quotient topologies visible without moving the patches.

The 3D view is intentionally lighter. It draws patch boundaries as translucent
surfaces, highlights both interface faces, places arrows along both native
tangential directions, and writes the full signed two-axis mapping beside the
plot. Solid and dashed face boundaries distinguish the minus and plus sides.
Shared and unshared topological vertices are colored and optionally labelled.
The result is intended for connectivity inspection rather than
publication-quality rendering of a complex volume decomposition.

Complete-vertex inspection and annotated plotting support 2D and 3D domains
whose physical dimension equals their logical dimension. The ordinary
geometric ``plot_domain`` paths remain available when ``topology=False``.

Selecting patch subdomains
==========================

``domain.get_subdomain(names)`` returns the topology induced by the requested
patch names.  An interface is retained when both of its sides belong to the
selection.  This includes repeated interfaces and interfaces from a patch to
itself.  When exactly one side is selected, that side becomes part of the new
exterior boundary.

For mapped patches, the returned domain receives the corresponding patch
mapping or ``MultiPatchMapping`` and a logical subdomain with matching
connectivity.  Interface objects and their names are retained, so orientation,
``axis_map``, ``InterfaceMapping``, and logical-interface identity are not
reconstructed or lost.

``Domain.join`` accepts atomic patch domains without existing connectivity.
Use ``get_subdomain`` to select from an existing joined domain; nested joined
domains are not implicitly flattened by another call to ``Domain.join``.

Exterior boundary construction
==============================

``Domain.join`` starts from every boundary face of every patch and removes all
faces used as an interface side.  The remaining faces form ``domain.boundary``.
Consequently a weak form integrated on ``domain.boundary`` sees only the
exterior boundary, not internal seams.

The join operation does not infer connectivity by comparing physical
coordinates.  Every interface must be declared.  Conversely, a declared face
is removed from the exterior boundary even if its two physical mappings do not
actually coincide.

What SymPDE validates
=====================

``Domain.join`` validates the structural contract:

* ``patches`` and ``connectivity`` have the expected container types;
* every input to ``Domain.join`` is an atomic, unconnected patch domain;
* patch objects, interior names, and serialized logical names are unique;
* all patch objects have a common logical dimension;
* patch indices are in range and patch objects belong to ``patches``;
* axes are valid and extremities are ``-1`` or ``+1``;
* every interface has an explicit dimension-appropriate orientation; and
* every orientation flag and sign belongs to ``{-1, +1}``.

Axes, extremities, and orientation flags must be integer values.  Floating
point values are rejected even when they are numerically equal to an integer,
and booleans are not interpreted as ``0`` or ``1``.

It does *not* establish geometric or discretization compatibility.  The caller
is responsible for ensuring that:

* the mapped physical faces coincide to the desired tolerance;
* patch maps are regular and have the intended physical dimension;
* the face parameterizations satisfy the declared orientation;
* a boundary face is not reused in a way that creates an unintended
  non-manifold topology; and
* later discrete trace meshes and spaces meet the requirements of the chosen
  Psydac operator.

This separation permits general quotient topologies, repeated patch-pair
interfaces, self-identifications, cross-axis gluing, reversed coordinates, and
3D face-axis permutations without baking a mesh policy into the symbolic
geometry.

Migration from the previous multipatch API
==========================================

The explicit orientation contract removes compatibility paths whose meaning
was ambiguous for cross-axis or three-dimensional interfaces.  Existing
multipatch callers should apply the following replacements:

.. list-table::
   :header-rows: 1

   * - Previous API or data
     - Current replacement
   * - ``Boundary.join(other, ornt=...)``
     - ``Boundary.join(other, orientation)``
   * - ``Interface(..., ornt=...)``
     - ``Interface(..., orientation)``
   * - ``interface.ornt``
     - ``interface.orientation``
   * - ``interface.axis``
     - ``interface.minus.axis`` and ``interface.plus.axis``
   * - ``Edge``
     - ``Boundary``
   * - interface dictionaries
     - ``(minus, plus, orientation[, name])`` tuples
   * - omitted interface orientation
     - an explicit dimension-appropriate orientation
   * - ``sympde.utilities.utils.plot_domain``
     - ``sympde.utilities.plot_domain``

The topology file-format migration is described in `Serialization`_.  In
particular, legacy interface entries without an explicit orientation are not
loaded implicitly.

Serialization
=============

``domain.todict()`` contains the interior, exterior boundary, and connectivity.
``domain.export('topology.h5')`` writes this dictionary as ``topology.yml`` in
an HDF5 file. The exterior boundary is always represented as a list: an empty
list for a fully closed domain, a one-element list for one exterior side, and a
longer list otherwise. The reader also accepts the historical single-side
dictionary representation.

Every serialized interface contains three entries: its minus side, plus side,
and explicit orientation. Its connectivity dictionary key is the interface
name and is preserved when the file is loaded. For example:

.. code-block:: yaml

   connectivity:
     A|B:
       - {patch: A, axis: '0', ext: '1', mapping: None, name: A_0}
       - {patch: B, axis: '1', ext: '-1', mapping: None, name: B_1}
       - [-1, 1, -1]

Load topology with ``Domain.from_file(filename)``.  Files written by the old
orientation-less format are deliberately rejected: assigning an implicit
orientation would be ambiguous and could silently change the topology.  Such
files must be regenerated or migrated. Files containing the experimental
``permutation``/``directions`` orientation mapping are likewise obsolete;
current files store the same compact values accepted by ``Domain.join``.
The loader derives the exterior boundary from the connectivity and rejects a
file when that result disagrees with boundary metadata that is present. Missing
or null legacy metadata is treated as unavailable and replaced by the derived
boundary.

SymPDE's ``Domain.export`` stores symbolic topology only.  A Psydac geometry
file additionally stores spline or NURBS knots, control points, weights, and
parallelizable mapping data; see Psydac's multipatch discretization tutorial.

Current scope
=============

The symbolic representation covers the following independently:

* 1D endpoint, 2D edge, and 3D face connections;
* unequal normal-axis indices on the two patches;
* all two 2D and all eight 3D face orientations;
* several interfaces between the same two patches;
* interfaces from a patch to itself; and
* arbitrary transitive corner incidence generated by those faces.

That is the scope of the *topology model*.  A downstream consumer may support
only a subset.  Psydac's standard variational and parallel interface paths
support signed cross-axis orientations for matching traces, while conforming
projectors and nested unequal trace couplings retain additional restrictions.
Those restrictions belong to Psydac and should not be encoded as assumptions
in a SymPDE domain.
