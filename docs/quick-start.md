# Quick-Start

- [Your first code using SymPDE](#sympde-poisson-0)
- [SymPDE concepts and their mathematical meaning](#sympde-concepts)
- [Examples](#sympde-examples)

<a id="sympde-poisson-0"></a>
## Your first code using SymPDE 
We first start by writing our first example using SymPDE. We consider the Poisson problem with homogeneous Dirichlet boundary conditions. 

```math
\begin{align}
  - \nabla^2 u = f \quad \text{in~$\Omega$}, \quad \quad 
  u = 0            \quad \text{on~$\partial \Omega$}. 
\end{align}
```
with $\Omega := (0,1)^2$ and

An $H^1$-conforming variational formulation of the previous problem reads
```math
\begin{align}
  \text{find $u \in V$ such that} \quad a(u,v) = l(v) \quad \forall v \in V,
\end{align}
```
where $V \subset H^1(\Omega)$, 
$a(u,v) := \int_{\Omega} \nabla u \cdot \nabla v ~ d\Omega$, and
$l(v) := \int_{\Omega} f v ~ d\Omega + \int_{\Gamma_N} g v ~ d\Gamma$.

The associated Python code can be found [here](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d/poisson_dir.py)

This very simple Python code reflects well the abstract mathematical framework needed for variational formulations.
The structure of the code is as follows,

1. Create a domain.
2. Create a space of *scalar* functions over the domain.
3. Create elements from this function space. These elements will denote the test and trial functions.
4. Create the Bilinear and Linear forms.
5. Create Essential Boundary Conditions.
6. Create the variational problem as an Equation.

<a id="sympde-concepts"></a>
## SymPDE concepts and their mathematical meaning

### Domain

sympde notation | Mathematical notion 
--- | ---
`Domain('Omega', dim=n)`  | an abstract domain $\Omega \subset \mathbb{R}^n$ 
`Domain('Omega', interiors, boundaries, connectivity)`  | an abstract domain $\Omega \subset \mathbb{R}^n$, represented by its interiors, boundaries and connectivity
`Interval('I', coordinate=x)`  | an abstract (interior) interval $I \subset \mathbb{R}$, with the coordinate $x$
`Union(A, B, C)`  | union of domains $A \cup B \cup C$
`ProductDomain(A, B)`  | product of two domains $A \times B$
`Indicator(A)`  | Indicator/characteristic function on $A$, *i.e.* $\mathbf{1}\_A$
`Line('Omega')`    | an abstract line segment, $\Omega \subset \mathbb{R}$, having one interval as interior domain
`Square('Omega')`  | an abstract square, $\Omega \subset \mathbb{R}^2$, having one interior domain as product of two intervals 
`Cube('Omega')`    | an abstract cube, $\Omega \subset \mathbb{R}^3$, having one interior domain as product of three intervals 
`NormalVector('n')` | Normal vector $\mathbf{n}$ 
`TangentVector('t')` | Tangent vector $\mathbf{t}$ 
`e = ElementDomain(Omega)` | an element of an abstract domain $e \in \mathcal{T}(\Omega)$ or $e = d\Omega$, where $\mathcal{T}(\Omega)$ is a given tessellation of $\Omega$ 
`DomainArea(Omega)` | Area of an abstract domain $\Omega$ 
`ElementArea(Omega)` | Area of an abstract element of a domain $\Omega$ 
`Area(A)` | Area of an expression of topological domain notions 

### Mapping

sympde notation | Mathematical notion 
--- | ---
`Mapping('F', n)` | a mapping functor $\mathbf{F}\_{\Omega} := \left( \mathbf{F}, \Omega \right): \Omega \rightarrow \mathbb{R}^n$
`F[i]`            | $i^{th}$ physical coordinate from $\{x,y,z\}$
`F.jacobian`      | the jacobian matrix $\mathcal{D}\_\mathbf{F}$
`F.det_jacobian`  | determinant of the jacobian matrix, $J_\mathbf{F}:= \mathrm{det} ~\mathcal{D}\_\mathbf{F}$
`F.covariant`     | the covariant matrix $\left( \mathcal{D}\_\mathbf{F} \right)^{-T}$
`F.contravariant` | the contravariant matrix $\frac{1}{J_\mathbf{F}} \mathcal{D}\_\mathbf{F}$ 
`F.hessian`       | the hessian matrix $\mathcal{H}\_\mathbf{F}$
`Covariant(F, v)` | action of the covariant matrix of $\mathbf{F}$ on $\mathbf{v}$, *i.e.* $\left( \mathcal{D}\_\mathbf{F} \right)^{-T} \mathbf{v}$
`Contravariant(F, v)` | action of the contravariant matrix of $\mathbf{F}$ on $\mathbf{v}$, *i.e.* $\frac{1}{J_\mathbf{F}} \mathcal{D}\_\mathbf{F} \mathbf{v}$

### Function spaces

sympde notation | Mathematical notion 
--- | ---
`ScalarFunctionSpace('V', Omega)`     | scalar function space $\mathcal{V}$
`VectorFunctionSpace('W', Omega)`     | vector function space $\mathbf{\mathcal{W}}$
`ProductSpace(V, W)}` or `V*W`        | product of spaces, *i.e.* $\mathcal{V} \times \mathcal{W}$
`element_of_space(V, 'v')`          | scalar function $v \in \mathcal{V}$ 
`element_of_space(W, 'w')`          | vector function $\mathbf{w} \in \mathcal{W}$  
`element_of_space(V*W, ['v', 'w'])` | $\left(v,\mathbf{w}\right) \in \mathcal{V} \times \mathcal{W}$ 

#### Function space types
sympde notation | Mathematical notion 
--- | ---
`H1SpaceType`    | $H^1$    
`HcurlSpaceType` | ${H}{\mbox{curl}}$ 
`HdivSpaceType`  | ${H}{\mbox{div}}$  
`L2SpaceType`    | $L^2$    

#### Typed function spaces

sympde notation | Mathematical notion 
--- | ---
`ScalarFunctionSpace('V0', Omega, kind='H1')`    | scalar function space $\mathcal{V}\_0 \subseteq H^1(\Omega)$    
`VectorFunctionSpace('V1', Omega, kind='Hcurl')` | vector function space $\mathcal{V}\_1 \subseteq H(\mbox{curl}, \Omega)$ 
`VectorFunctionSpace('V2', Omega, kind='Hdiv') ` | vector function space $\mathcal{V}\_2 \subseteq H(\mbox{div}, \Omega)$  
`ScalarFunctionSpace('V3', Omega, kind='L2')`    | scalar function space $\mathcal{V}\_3 \subseteq L^2(\Omega)$    

### Projections

sympde notation | Mathematical notion 
--- | ---
`P_V := Projector(V)`                 | Projector onto the function space $\mathcal{V}$  
`P_V(expr)`                           | Projection of the expression `expr` onto the function space $\mathcal{V}$  
`Pi_V := Projector(V, 'commuting')`   | Commuting projector onto the typed function space $\mathcal{V}$  
`Pi_V(expr)`                          | Commuting projection of the expression `expr` onto the typed function space $\mathcal{V}$  

### Atomic variables

sympde notation | Mathematical notion 
--- | ---
`ScalarFunction(V, 'u')`    |  $u \in \mathcal{V}$
`VectorFunction(W, 'w')`    |  $\mathbf{w} \in \mathcal{W}$
`Constant('mu', real=True)` | constant number $\mu \in \mathbb{R}$ 
`int`                       | integer number 
`float`                     | floating-point number

### Algebraic operators 

sympde notation | Mathematical notion 
--- | ---
`dot(u, v)`    | $\mathbf{u} \cdot \mathbf{v}$   
`inner(u, v)`  | $\mathbf{u} : \mathbf{v}$       
`cross(u, v)`  | $\mathbf{u} \times \mathbf{v}$  
`outer(u, v)`  | $\mathbf{u} \otimes \mathbf{v}$  

### Differential operators 

boldface font is used for vector functions/expressions

sympde notation | Mathematical notion 
--- | ---
`dx(u)`, `dy(u)` or `dz(u)`   |  $\partial_x u$, $\partial_y u$ or $\partial_z u$                               
`grad(u)`                     |  $\nabla u$ or $\nabla \mathbf{u}$                                               
`div(u)`                      |  $\nabla \cdot \mathbf{u}$                                                         
`curl(u)`                     |  $\nabla \times \mathbf{u}$                                                        
`rot(u) `                     |  2D rotational $\nabla \times u$                                                   
`convect(a, u)`               |  $\left( \mathbf{a} \cdot \nabla \right) \mathbf{u}$                               
`trace(u)`                    |  trace $\gamma(u)$                                                                 
`Dn(u)`                       |  normal derivative $\partial_{\mathbf{n}} u$                                       
`D(u)`                        |  Strain tensor $\frac{1}{2}\left(\nabla \mathbf{u} + \nabla \mathbf{u}^T \right) $ 
`laplace(u)`                  |  Laplace $\Delta u$                                                                
`hessian(u)`                  |  Hessian $H(u)$                                                                    
`bracket(u,v)`                |  Poisson Bracket $[u,v]$                                                           
              
### Additional operators

sympde notation | Mathematical notion 
--- | ---
`jump(u)`      |  jump of $u$, *i.e.* $[u]$                            
`avg(u)`       |  average of $u$, *i.e.* $\langle u \rangle$           
`conv(K,u)`    |  convolution of $u$ with a kernel $K$, *i.e.* $K * u$ 

### Integral operators 

sympde notation | Mathematical notion 
--- | ---
`DomainIntegral(f)`   | integral over a domain, *i.e.* $(f, \Omega) \mapsto \int_{\Omega} f ~d\Omega$
`BoundaryIntegral(f)` | integral over a boundary, *i.e.* $(f, \Gamma) \mapsto \int_{\Gamma} f ~d\partial\Omega$
`PathIntegral(F)`     | integral over an oriented path, *i.e.* $(\mathbf{F}, C) \mapsto \int_{C} \mathbf{F} \cdot d\mathbf{s}$
`SurfaceIntegral(F)`  | integral over an oriented surface, *i.e.* $(\mathbf{F}, S) \mapsto \int_{S} \mathbf{F} \cdot d\mathbf{S}$

<a id="sympde-examples"></a>
## Examples

