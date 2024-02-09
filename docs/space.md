# Function Space concepts

- [BasicFunctionSpace](#BasicFunctionSpace)
- [ScalarFunctionSpace](#ScalarFunctionSpace)
- [VectorFunctionSpace](#VectorFunctionSpace)
- [ProductSpace](#ProductSpace)
- [Functions](#Functions)
- [Derham](#Derham)
- [Trace](#Trace)

SymPDE provides two Python classes to describe scalar and vector function spaces, respectively. There is no notion of a discrete representation; for example, we do not need to mention that a space is a Brezzi-Douglas-Marini (BDM) space. In fact, UFL function spaces can be seen as parametric types having as a basic type FunctionSpace or VectorFunctionSpace. SymPDE only needs to know if an element of the space can be indexed or not. A BDM space would then be identified by an annotation added by a third party library to uniquely define a function space at the discrete level.

<a id="BasicFunctionSpace"></a>
## BasicFunctionSpace 

<a id="ScalarFunctionSpace"></a>
## ScalarFunctionSpace 

<a id="VectorFunctionSpace"></a>
## VectorFunctionSpace 

<a id="ProductSpace"></a>
## ProductSpace 

<a id="Derham"></a>
## Derham 

<a id="Functions"></a>
## Functions

### ScalarFunction
### IndexedVectorFunction
### VectorFunction

<a id="Trace"></a>
## Trace 

### trace_0
### trace_1

## Utilities
### element_of and elements_of

## Projection (Experimental)

## Projector (Experimental) 

 
