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
`F.det\_jacobian` | determinant of the jacobian matrix, $J_\mathbf{F}:= \mathrm{det} ~\mathcal{D}\_\mathbf{F}$
`F.covariant`     | the covariant matrix $\left( \mathcal{D}\_\mathbf{F} \right)^{-T}$
`F.contravariant` | the contravariant matrix $\frac{1}{J_\mathbf{F}} \mathcal{D}\_\mathbf{F}$ 
`F.hessian`       | the hessian matrix $\mathcal{H}\_\mathbf{F}$
`Covariant(F, v)` | action of the covariant matrix of $\mathbf{F}$ on $\mathbf{v}$, \textit{i.e.} $\left( \mathcal{D}\_\mathbf{F} \right)^{-T} \mathbf{v}$
`Contravariant(F, v)` | action of the contravariant matrix of $\mathbf{F}$ on $\mathbf{v}$, \textit{i.e.} $\frac{1}{J_\mathbf{F}} \mathcal{D}\_\mathbf{F} \mathbf{v}$
