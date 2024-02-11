# Function Space concepts

- [BasicFunctionSpace](#BasicFunctionSpace)
- [ScalarFunctionSpace](#ScalarFunctionSpace)
- [VectorFunctionSpace](#VectorFunctionSpace)
- [Typed Spaces](#typed-spaces)
- [ProductSpace](#ProductSpace)
- [Functions](#Functions)
- [Derham](#Derham)
- [Trace](#Trace)

SymPDE provides two Python classes to describe scalar and vector function spaces, respectively. There is no notion of a discrete representation; for example, we do not need to mention that a space is a Brezzi-Douglas-Marini (BDM) space. In fact, these kind of function spaces can be seen as parametric types having as a basic type **ScalarFunctionSpace** or **VectorFunctionSpace**. SymPDE only needs to know if an element of the space can be indexed or not. A BDM space would then be identified by an annotation added by a third party library to uniquely define a function space at the discrete level.

<a id="BasicFunctionSpace"></a>
## BasicFunctionSpace 

This represents the base class of all our function spaces. It does not reflect a continuous space, but it is more a **type** at the formal/abstract level.

<a id="ScalarFunctionSpace"></a>
## ScalarFunctionSpace 

```python
>>> domain = Domain('Omega', dim=2)
>>> V1 = ScalarFunctionSpace('V1', domain)
```

<a id="VectorFunctionSpace"></a>
## VectorFunctionSpace 

```python
>>> domain = Domain('Omega', dim=2)
>>> W1 = VectorFunctionSpace('W1', domain)
```

<a id="typed-spaces"></a>
## Typed Spaces
```python
>>> domain = Domain('Omega', dim=2)
>>> H1    = ScalarFunctionSpace('V0', domain, kind='H1')
>>> Hcurl = VectorFunctionSpace('V1', domain, kind='Hcurl')
>>> Hdiv  = VectorFunctionSpace('V2', domain, kind='Hdiv')
>>> L2    = ScalarFunctionSpace('V3', domain, kind='L2')
>>> assert(H1.kind    == H1Space)                
>>> assert(Hcurl.kind == HcurlSpace)
>>> assert(Hdiv.kind  == HdivSpace)
>>> assert(L2.kind    == L2Space)
>>> assert(V.kind     == UndefinedSpace)
>>> assert(W.kind     == UndefinedSpace)
>>> assert(H1.regularity    > L2.regularity)
>>> assert(H1.regularity    > Hcurl.regularity)
>>> assert(Hcurl.regularity > L2.regularity)
```

<a id="ProductSpace"></a>
## ProductSpace 

```python
>>> domain = Domain('Omega', dim=2)
>>> V1 = ScalarFunctionSpace('V1', domain)
>>> V2 = ScalarFunctionSpace('V2', domain)
>>> W1 = VectorFunctionSpace('W1', domain)
>>> X = V1 * V2 * W1
>>> # or 
>>> X = ProductSpace(V1, V2, W1)
```

<a id="Derham"></a>
## Derham 

<a id="Functions"></a>
## Functions

### ScalarFunction
A scalar function is necessary an element of the **ScalarFunctionSpace**.

### IndexedVectorFunction
### VectorFunction
A vector function is necessary an element of the **VectorFunctionSpace**. It is an extension of the **IndexedVectorFunction**, which means that one can access to its components.

<a id="Trace"></a>
## Trace 

### trace_0
### trace_1

## Utilities
### element_of and elements_of

## Projection (Experimental)

## Projector (Experimental) 

 
