# Topology concepts


- [Domain](#domain)
- [PeriodicDomain](#periodic-domain)
- [NCube](#ncube-domain)
- [InteriorDomain](#interior-domain)
- [NCubeInterior](#ncube-interior-domain)
- [Additional geometry objects](#geometry-addons)

<a id="domain"></a>
## Domain

Mathematicaly speaking, a domain $\Omega$ is defined as an interior $\a{\Omega}$ and a boundary $\partial \Omega$.
A domain is defined by its **interior** and its **boundaries**. For instance, a domain without a boundary is periodic, since, **sympde** does not consider infinite domains.

In addition to these two attributes, a **Domain** can also be defined as the image of a *logical domain* using a *mapping*.

Finally, a **Domain** can have subdomains, in which case, a **connectivity** must be provided as well.

The simplest **Domain** is the **NCube** object. For instance, **Line**, **Square** and **Cube** are specific cases of **NCube**.

Let's take a look at the **Line** object.

```python
>>> from sympde import topology
>>> domain = topology.Line()
>>> print(type(domain))
sympde.topology.domain.Line
```

The parent class of **Line** is **NCube** as we can see in our example

```python
>>> print(type(topology.Line.__bases__))
sympde.topology.domain.NCube
```

In fact, the objects **Line**, **Square** and **Cube** are subclasses of **NCube** which is a subclasse of **Domain** that is an extension of the **BasicDomain** object.


### Properties

- **name**: returns the domain's name
- **interior**: returns the interior of the domain
- **boundary**: returns the boundary of the domain
- **mapping**: returns the mapping, if associated to the domain
- **logical_domain**: returns the logical domain if a mapping is used
- **connectivity**: returns information about the interfaces of the subdomains
- **dim**: dimension of the space
- **dtype**: a dictionary that contains information about the domain
- **interfaces**: returns the interfaces of the subdomains 
- **corners**: returns the corners of the domain

**TODO** add inherited properties

### Examples

#### Example 1

```python
>>> from sympde.topology import InteriorDomain, Boundary, Domain

>>> Omega_1 = InteriorDomain('Omega_1', dim=2)
>>> Omega_2 = InteriorDomain('Omega_2', dim=2)

>>> Gamma_1 = Boundary('Gamma_1', Omega_1)
>>> Gamma_2 = Boundary('Gamma_2', Omega_2)
>>> Gamma_3 = Boundary('Gamma_3', Omega_2)

>>> Omega = Domain('Omega',
                   interiors=[Omega_1, Omega_2],
                   boundaries=[Gamma_1, Gamma_2, Gamma_3])

>>> print( Omega.dim )
2
>>> print( len(Omega.interior) )
2
>>> print( len(Omega.boundary) )
3
```

#### Example 2

```python
>>> from sympde.topology import InteriorDomain, Boundary, Domain

>>> Omega_1 = InteriorDomain('Omega_1', dim=2)

>>> Gamma_1 = Boundary('Gamma_1', Omega_1)
>>> Gamma_2 = Boundary('Gamma_2', Omega_1)

>>> Omega = Domain('Omega',
                   interiors=[Omega_1],
                   boundaries=[Gamma_1, Gamma_2])

>>> assert(Omega.boundary == Union(Gamma_1, Gamma_2))
>>> assert(Omega.boundary.complement(Gamma_1) == Gamma_2)
>>> assert(Omega.boundary - Gamma_1 == Gamma_2)
```

<a id="periodic-domain"></a>
## PeriodicDomain

Similar to the **Domain** object, a **PeriodicDomain** is defined by giving a **Domain** and the associated **periods**.

### Properties

- **domain**: returns the domain
- **periods**: returns the flags for each axis
- **boundary**: return the boundary of the domain
- **dim** returns the dimension of the domain
- **coordinates** return the coordinates symbols

<a id="ncube-domain"></a>
## NCube

An **NCube** is a **Domain** that models hypercubes in $\mathbb{R}^n$. It is defined by specifying its **name**, **dim**, **min_coords** and **max_coords**.

The interior of a **NCube** is defined as **NCubeInterior**.
To make things easier, we provide specific classes in the case where $n \in \{ 1, 2, 3 \}$.

### Line
This corresponds to the case where $n=1$, meaning it is $(x_{min}, x_{max})$.

#### Properties
- **bounds**: the bounds for the line, i.e. $x_{min}$ and $x_{max}$.

### Square
This corresponds to the case where $n=2$, meaning it is $(x_{min}, x_{max}) \times (y_{min}, y_{max})$.

#### Properties
- **bounds1**: the bounds in the first axis, i.e. $x_{min}$ and $x_{max}$.
- **bounds2**: the bounds in the first axis, i.e. $y_{min}$ and $y_{max}$.

### Cube
This corresponds to the case where $n=3$, meaning it is $(x_{min}, x_{max}) \times (y_{min}, y_{max}) \times (z_{min}, z_{max})$.

#### Properties
- **bounds1**: the bounds in the first axis, i.e. $x_{min}$ and $x_{max}$.
- **bounds2**: the bounds in the first axis, i.e. $y_{min}$ and $y_{max}$.
- **bounds3**: the bounds in the first axis, i.e. $z_{min}$ and $z_{max}$.


<a id="interior-domain"></a>
## InteriorDomain

**TODO**

### Examples
#### Example 1

In this example, we create two interior domains and then we build their union.

```python
>>> from sympde.topology import InteriorDomain, Union

>>> D1 = InteriorDomain('D1', dim=2)
>>> D2 = InteriorDomain('D2', dim=2)

>>> # Union is commutative
>>> assert( Union(D2, D1) == Union(D1, D2) )

>>> D = Union(D1, D2)
>>> # the union will be of dimension 2, and has a length of 2
>>> print(D.dim)
2
>>> print(len(D))
2
>>> # we see it clearly when converting it to a dictionary
>>> print(D.todict())
[{'name': 'D1', 'mapping': 'None'}, {'name': 'D2', 'mapping': 'None'}]
```


<a id="ncube-interior-domain"></a>
## NCubeInterior

This class defines the interior of a domain which is **NCube**.

### Properties

- **min_coords**
- **max_coords**
- **boundary**


<a id="geometry-addons"></a>
## Additional geometry objects

### Vectors

#### BoundaryVector

is an extension of the **IndexedBase** object from **sympy**.

#### NormalVector

defines the normal vector to boundaries. It is an extension of the **BoundaryVector** object.

#### MinusNormalVector

defines the inward normal vector to a boundary. It is an extension of the **NormalVector** object.

#### MinusNormalVector

defines the outward normal vector to a boundary. It is an extension of the **NormalVector** object.

#### TangentVector

defines the tangent vector to boundaries. It is an extension of the **BoundaryVector** object.

### ElementDomain

defines an element of the domain.

### BasicArea

A basic class defining the concept of an area. It is an extension of the **AtomicExpr** object from **sympy**, this way we can use it in arithmetic expressions.

### DomainArea

defines the area of a **Domain**.

### ElementArea

defines the area of an element of a **Domain**.

### Geometry Operations

**TODO**

#### BasicGeometryOperator

**TODO**

#### Area

**TODO**
