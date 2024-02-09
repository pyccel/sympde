# Topological concepts

## Domain concepts and their mathematical meaning

sympde notation | Mathematical notion 
--- | ---
`Domain('Omega', dim=n)`  | an abstract domain $\Omega \subset \mathbb{R}^n$ 
`Domain('Omega', interiors, boundaries, connectivity)`  | an abstract domain $\Omega \subset \mathbb{R}^n$, represented by its interiors, boundaries and connectivity
`Interval('I', coordinate=x)`  | an abstract (interior) interval $I \subset \mathbb{R}$, with the coordinate $x$
`Union(A, B, C)`  | union of domains $A \cup B \cup C$
`ProductDomain(A, B)`  | product of two domains $A \times B$
`Indicator(A)`  | Indicator/characteristic function on $A$, \textit{i.e.} $\mathbf{1}\_A$
`Line('Omega')`    | an abstract line segment, $\Omega \subset \mathbb{R}$, having one interval as interior domain
`Square('Omega')`  | an abstract square, $\Omega \subset \mathbb{R}^2$, having one interior domain as product of two intervals 
`Cube('Omega')`    | an abstract cube, $\Omega \subset \mathbb{R}^3$, having one interior domain as product of three intervals 
`NormalVector('n')` | Normal vector $\mathbf{n}$ 
`TangentVector('t')` | Tangent vector $\mathbf{t}$ 
`e = ElementDomain(Omega)` | an element of an abstract domain $e \in \mathcal{T}(\Omega)$ or $e = d\Omega$, where $\mathcal{T}(\Omega)$ is a given tessellation of $\Omega$ 
`DomainArea(Omega)` | Area of an abstract domain $\Omega$ 
`ElementArea(Omega)` | Area of an abstract element of a domain $\Omega$ 
`Area(A)` | Area of an expression of topological domain notions 

## Mapping concepts their mathematical meaning  

sympde notation | Mathematical notion 
--- | ---
`Mapping('F', n)` | a mapping functor $\mathbf{F}\_{\Omega} := \left( \mathbf{F}, \Omega \right): \Omega \rightarrow \mathbb{R}^n$
`F[i]`            | $i^{th}$ physical coordinate from $\{x,y,z\}$
`F.jacobian`      | the jacobian matrix $\mathcal{D}\_\mathbf{F}$
`F.det_jacobian`  | determinant of the jacobian matrix, $J_\mathbf{F}:= \mathrm{det} ~\mathcal{D}\_\mathbf{F}$
`F.covariant`     | the covariant matrix $\left( \mathcal{D}\_\mathbf{F} \right)^{-T}$
`F.contravariant` | the contravariant matrix $\frac{1}{J_\mathbf{F}} \mathcal{D}\_\mathbf{F}$ 
`F.hessian`       | the hessian matrix $\mathcal{H}\_\mathbf{F}$
`Covariant(F, v)` | action of the covariant matrix of $\mathbf{F}$ on $\mathbf{v}$, \textit{i.e.} $\left( \mathcal{D}\_\mathbf{F} \right)^{-T} \mathbf{v}$
`Contravariant(F, v)` | action of the contravariant matrix of $\mathbf{F}$ on $\mathbf{v}$, \textit{i.e.} $\frac{1}{J_\mathbf{F}} \mathcal{D}\_\mathbf{F} \mathbf{v}$

## Function spaces and their mathematical meaning

sympde notation | Mathematical notion 
--- | ---
`ScalarFunctionSpace('V', Omega)`     | scalar function space $\mathcal{V}$
`VectorFunctionSpace('W', Omega)`     | vector function space $\mathbf{\mathcal{W}}$
`ProductSpace(V, W)}` or `V*W`        | product of spaces, \textit{i.e.} $\mathcal{V} \times \mathcal{W}$
`element\_of\_space(V, 'v')`          | scalar function $v \in \mathcal{V}$ 
`element\_of\_space(W, 'w')`          | vector function $\mathbf{w} \in \mathcal{W}$  
`element\_of\_space(V*W, ['v', 'w'])` | $\left(v,\mathbf{w}\right) \in \mathcal{V} \times \mathcal{W}$ 

### Function space types
sympde notation | Mathematical notion 
--- | ---
`H1SpaceType`    | $H^1$    
`HcurlSpaceType` | ${H}{\mbox{curl}}$ 
`HdivSpaceType`  | ${H}{\mbox{div}}$  
`L2SpaceType`    | $L^2$    

### Declaration of typed function spaces

sympde notation | Mathematical notion 
--- | ---
`ScalarFunctionSpace('V0', Omega, kind='H1')`    | scalar function space $\mathcal{V}\_0 \subseteq H^1(\Omega)$    
`VectorFunctionSpace('V1', Omega, kind='Hcurl')` | vector function space $\mathcal{V}\_1 \subseteq H(\mbox{curl}, \Omega)$ 
`VectorFunctionSpace('V2', Omega, kind='Hdiv') ` | vector function space $\mathcal{V}\_2 \subseteq H(\mbox{div}, \Omega)$  
`ScalarFunctionSpace('V3', Omega, kind='L2')`    | scalar function space $\mathcal{V}\_3 \subseteq L^2(\Omega)$    

## Projections and their mathematical meaning

sympde notation | Mathematical notion 
--- | ---
`P_V := Projector(V)`                 | Projector onto the function space $\mathcal{V}$  
`P_V(expr)`                           | Projection of the expression \texttt{expr} onto the function space $\mathcal{V}$  
`Pi_V := Projector(V, 'commuting')`   | Commuting projector onto the typed function space $\mathcal{V}$  
`Pi_V(expr)`                          | Commuting projection of the expression \texttt{expr} onto the typed function space $\mathcal{V}$  


