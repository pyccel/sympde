# Examples

- [Linear Problems](#linear-problems)
  - [The Poisson problem](#linear-poisson)
  - [The Biharmonic problem](#linear-biharmonic)
  - [The Vector Poisson problem](#linear-vector-poisson)
  - [The advection-diffusion equation](#linear-advection-diffusion)
  - [The Stabilized advection-diffusion equation](#linear-stabilized-advection-diffusion)
  - [Mixed FEM for the Poisson problem](#linear-mixed-poisson)
  - [Mixed FEM for the Stokes problem](#linear-mixed-stokes)
  - [Elliptic-curl Problem](#linear-elliptic-curl)
  - [Elliptic-div Problem](#linear-elliptic-div)
  - [Linear Elasticity Problem](#linear-elasticity)
- [Nonlinear Problems](#nonlinear-problems)
  - [Nonlinear 2D Poisson Problem](#nonlinear-poisson-2d)
  - [1D Burgers Problem](#nonlinear-burgers-1d)
  - [The streamfunction-velocity formulation of the steady-state Navier-Stokes equations for incompressible fluids](#nonlinear-navier-stokes-streamfunction-velocity)
- [Subdomains](#subdomains)
  - [The Poisson problem using Nitsche method](#subdomains-linear-poisson)


<a id="linear-problems"></a>
## Linear Problems

<a id="linear-poisson"></a>
### The Poisson problem

As a first example, we consider the Poisson equation

$$
\begin{align}
  \- \nabla^2 u = f \quad \text{in~$\Omega$}, \quad \quad 
  u = 0            \quad \text{on~$\Gamma_D$}, \quad \quad 
  \partial_n u = g \quad \text{on~$\Gamma_N$}.
\end{align}
$$

An $H^1$-conforming variational formulation of reads

$$
\begin{align}
  \text{find $u \in V$ such that} \quad a(u,v) = l(v) \quad \forall v \in V,
\end{align}
$$

where 

- $V \subset H^1(\Omega)$, 
- $a(u,v) := \int_{\Omega} \nabla u \cdot \nabla v ~ d\Omega$,
- $l(v) := \int_{\Omega} f v ~ d\Omega + \int_{\Gamma_N} g v ~ d\Gamma$.

#### The Poisson equation with weak imposition of Dirichlet conditions

In this example, we consider the following Poisson problem

$$
\begin{align}
  \- \nabla^2 u = f \quad \text{in~$\Omega$}, \quad \quad 
  u = g            \quad \text{on~$\partial \Omega$}.
\end{align}
$$

The variational formulation writes

$$
\begin{align}
  \text{find $u \in V$ such that} \quad a(u,v) = l(v) \quad \forall v \in V,
\end{align}
$$

where 

- $V \subset H^1(\Omega)$, 
- $a(u,v) := \int_{\Omega} \nabla u \cdot \nabla v ~ d\Omega + \int_{\partial \Omega} \left( \kappa u v - u \partial_n v - v \partial_n u \right)~d\partial\Omega $,
- $l(v) := \int_{\Omega} f v ~ d\Omega + \int_{\partial\Omega} \left( \kappa g v - g \partial_n v \right) ~ d\Gamma$.


#### Examples

- [Poisson on a square domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_poisson_dir0.py)
- [Poisson on a square domain with Homogeneous and Neumann Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_poisson_dir0_neu.py)
- [Poisson on a square domain with Homogeneous, Inhomogeneous and Neumann Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_poisson_dir_neu.py)
- [Poisson equation on a square with weak imposition of Dirichlet conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_poisson_nitsche_diri.py)

<a id="linear-biharmonic"></a>
### The Biharmonic problem

We consider the (inhomogeneous) biharmonic equation with homogeneous essential boundary conditions,

$$
\begin{align}
  \nabla^2 \nabla^2  u = f       \quad \text{in $\Omega$}         , \quad \quad 
  u = 0                          \quad \text{on $\partial \Omega$}, \quad \quad
  \nabla u \cdot \mathbf{n} = 0  \quad \text{on $\partial \Omega$}.
\end{align}
$$

An $H^2$-conforming variational formulation reads 

$$
\begin{align}
  \text{find $u \in V$ such that} \quad a(u,v) = l(v) \quad \forall v \in V,
\end{align}
$$

where 

- $V \subset H^2(\Omega)$, 
- $a(u,v) := \int_{\Omega} \nabla^2 u ~ \nabla^2 v ~ d\Omega$,
- $l(v) := \int_{\Omega} f v ~ d\Omega$.

#### Examples

- [Biharmonic equation on a square domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_biharmonic_dir0.py)

<a id="linear-vector-poisson"></a>
### Vector Poisson's equation 

In this example we consider the vector Poisson equation with homogeneous Dirichlet boundary conditions:

$$
\begin{align}
  \- \nabla^2 \mathbf{u} = \mathbf{f} \quad \mbox{in} ~ \Omega, \quad \quad 
  \mathbf{u} = 0            \quad \mbox{on} ~ \partial \Omega.
\end{align}
$$

The corresponding variational formulation, using $\mathbf{H}^1$ formulation, *i.e.* all components are in $H^1$, reads 

$$
\begin{align}
  \text{find $\mathbf{u} \in V$ such that} \quad 
  a(\mathbf{u},\mathbf{v}) = l(\mathbf{v}) \quad \forall \mathbf{v} \in V,
\end{align}
$$

where 

- $V \subset \mathbf{H}\_0^1(\Omega)$, 
- $a(\mathbf{u},\mathbf{v}) := \int_{\Omega} \nabla \mathbf{u} : \nabla \mathbf{v} ~ d\Omega$,
- $l(\mathbf{v}) := \int_{\Omega} \mathbf{f} \cdot \mathbf{v} ~ d\Omega$.

#### Examples

- [Vector Poisson on a cube domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/3d_vector_poisson_dir0.py)

<a id="linear-advection-diffusion"></a>
### Advection-diffusion equation

We consider the advection-diffusion problem consisting of finding a scalar-valued function $u$ such that

$$
\begin{align}
  \begin{cases}
    - \kappa \nabla^2 u + \mathbf{b} \cdot \nabla u = f &\text{in $\Omega$}, \\
    u = 0                                               &\text{on $\partial \Omega$}.
  \end{cases}
\end{align}
$$

A variational formulation reads

$$
\begin{align}
  \text{find $u \in V$ such that} \quad a(u,v) = l(v) \quad \forall v \in V,
\end{align}
$$

where 

- $V \subset H_0^1(\Omega)$, 
- $a(u,v) := \int_{\Omega} \left( \kappa \nabla u \cdot \nabla v + \left( \mathbf{b} \cdot \nabla u \right) v \right) d\Omega$,
- $l(v) := \int_{\Omega} f v~d\Omega$.

#### Examples

- [Advection Diffusion problem on a square domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_advection_diffusion_dir0.py)

<a id="linear-stabilized-advection-diffusion"></a>
### Stabilized advection-diffusion equation

In the sequel we consider two different stabilization methods, leading to formulations of the form

$$
\begin{align}
  \text{find $u \in V$ such that} \quad 
  a(u,v) + a_s(u,v) = l(v) + l_s(v) \quad \forall v \in V,
\end{align}
$$

where $a_s$ and $l_s$ are additional terms depending on the stabilization method.

The streamline upwind Petrov-Galerkin (SUPG) method is given by

$$
\begin{align}
    a_s(u,v) = \int_{\Omega} \tau ~ L(u) ~ \left( \mathbf{b} \cdot \nabla v \right) ~d\Omega 
    \quad \mbox{and} \quad
    l_s(v) = \int_{\Omega} \tau ~ f ~ \left( \mathbf{b} \cdot \nabla v \right) ~d\Omega 
  && \text{[SUPG]}
\end{align}
$$

while the Galerkin least squares (GLS) method is as follows, 

$$
\begin{align}
    a_s(u,v) = \int_{\Omega} \tau ~ L(u) ~ L(v) ~d\Omega 
    \quad \mbox{and} \quad
    l_s(v) = \int_{\Omega} \tau ~ f ~ L(v) ~d\Omega 
  && \text{[GLS]}
\end{align}
$$

where we introduced the differential operator $L(u) := - \kappa \nabla^2 u + \mathbf{b} \cdot \nabla u$.

Notice that in these formulations $\tau$ stands for a piecewise function that is constant on each element of the tessellation associated to the computational domain $\Omega$.

#### Examples

- [Stabilized Advection Diffusion problem on a square domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_advection_diffusion_stabilized_dir0.py)


<a id="linear-mixed-poisson"></a>
### Mixed FEM for the Poisson problem

Let $\Omega \subset \mathbb{R}^3$ and consider the Poisson problem

$$
\begin{align}
  \left\{ 
  \begin{array}{clr}
    -\Delta p & =f & ,~\Omega    \\
    p         & =0 & ,~\partial \Omega
  \end{array} \right.
\end{align}
$$

Using that $\Delta p = \nabla\cdot\nabla p$, we set $ \mathbf{u}=\nabla p$, then the Poisson equation can be written equivalently

$$ \mathbf{u}=-\nabla p, ~~~ \nabla\cdot \mathbf{u}= f.$$

#### First mixed formulation of the Poisson problem

Instead of having one unknown, we now have two, along with the above two equations.
In order to get a mixed variational formulation, we first take the dot product of the first one by $ \mathbf{v}$ and integrate by parts

$$
\begin{align}
\int_{\Omega} \mathbf{u}\cdot \mathbf{v}~\mathrm{d} \mathbf{x} -\int_{\Omega} p ~ \nabla\cdot \mathbf{v}~\mathrm{d} \mathbf{x} + \int_{\partial\Omega} p ~  \mathbf{v}\cdot \mathbf{n}~\mathrm{d} \sigma = \int_{\Omega} \mathbf{u}\cdot \mathbf{v}~\mathrm{d} \mathbf{x} -\int_{\Omega} p ~ \nabla\cdot \mathbf{v}~\mathrm{d} \mathbf{x}=0,
\end{align}
$$

using $p=0$ as a natural boundary condition. Then multiplying the second equation by $q$ and integrating yields

$$
\begin{align}
\int_{\Omega} \nabla\cdot\mathbf{u} ~ q ~\mathrm{d} \mathbf{x} = \int_{\Omega} f q ~\mathrm{d} \mathbf{x}.
\end{align}
$$

No integration by parts is necessary here. And we thus get the following mixed variational formulation:

Find $(\mathbf{u},p) \in H(\mathrm{div},\Omega)\times L^2(\Omega)$ such that

$$
\begin{align}
\left\{ 
\begin{array}{llll}
  \int_{\Omega} \mathbf{u}\cdot \mathbf{v}~\mathrm{d} \mathbf{x} &- \int_{\Omega} p ~ \nabla\cdot \mathbf{v}~\mathrm{d} \mathbf{x} &=0, & \forall \mathbf{v}\in H(\mathrm{div},\Omega) \\
  \int_{\Omega} \nabla\cdot\mathbf{u} ~ q ~\mathrm{d} \mathbf{x} &  &= \int_{\Omega} f q ~\mathrm{d} \mathbf{x}, & \forall q\in L^2(\Omega)
\end{array} \right.
\end{align}
$$

#### Examples

- [Poisson problem with mixed (first) formulation on a square domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_poisson_mixed_v1_dir0.py)

#### Second mixed formulation of the Poisson problem

Here, we get an alternative formulation by not integrating by parts, the mixed term in the first formulation but in the second. The first formulation simply becomes

$$
\begin{align}
\int_{\Omega} \mathbf{u}\cdot \mathbf{v}~\mathrm{d} \mathbf{x} +\int_{\Omega} \nabla p \cdot \mathbf{v}~\mathrm{d} \mathbf{x}=0,
\end{align}
$$

and the second, removing immediately the boundary term due to the essential boundary condition $q=0$

$$
\begin{align}
\int_{\Omega}\nabla \cdot\mathbf{u}  ~ q ~\mathrm{d} \mathbf{x} = -\int_{\Omega}  \mathbf{u} \cdot \nabla q  ~\mathrm{d} \mathbf{x} = \int_{\Omega} f q ~\mathrm{d} \mathbf{x},
\end{align}
$$

which leads to the variational formulation

Find $(\mathbf{u},p) \in L^2(\Omega)^3 \times H^1_0(\Omega)$ such that

$$
\begin{align}
\left\{ 
\begin{array}{llll}
  \int_{\Omega} \mathbf{u}\cdot \mathbf{v}~\mathrm{d} \mathbf{x} &+ \int_{\Omega} \nabla p \cdot \mathbf{v}~\mathrm{d} \mathbf{x} &=0, & \forall \mathbf{v}\in L^2(\Omega)^3 \\
  \int_{\Omega}  \mathbf{u} \cdot \nabla q  ~\mathrm{d} \mathbf{x} & & = -\int_{\Omega} f q ~\mathrm{d} \mathbf{x}, & \forall q\in H^1_0(\Omega)
\end{array} \right.
\end{align}
$$


Note that this formulation actually contains the classical variational formulation for the Poisson equation. Indeed for $q\in H^1_0(\Omega)$, $\nabla q \in L^2(\Omega)^3$ can be used as a test function in the first equation. And plugging this into the second we get

$$
\int_{\Omega}  \nabla p \cdot \nabla q  ~\mathrm{d} \mathbf{x}  = \int_{\Omega} f q ~\mathrm{d} \mathbf{x}, \quad \forall q\in H^1_0(\Omega).
$$

#### Examples

- [Poisson problem with mixed (second) formulation on a square domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_poisson_mixed_v2_dir0.py)

<a id="linear-mixed-stokes"></a>
### Mixed FEM for the Stokes problem 

We consider now the Stokes problem for the steady-state modelling of an incompressible fluid

$$
\begin{align}
  \left\{
    \begin{array}{rl}
      - \nabla^2 \mathbf{u} + \nabla p = \mathbf{f} & \mbox{in} ~ \Omega,  \\
      \nabla \cdot \mathbf{u}          = 0   & \mbox{in} ~ \Omega,  \\
      \mathbf{u}               = 0   & \mbox{on} ~ \partial \Omega,
    \end{array}
    \right.
\end{align}
$$

#### First mixed formulation of the Stokes problem

For the variational formulation, we take the dot product of the first equation with $v$ and integrate over the whole domain

$$
\int_{\Omega} (-\Delta \mathbf{u} + \nabla p)\cdot \mathbf{v} ~\mathrm{d} \mathbf{x} = \int_{\Omega}\nabla \mathbf{u}:\nabla \mathbf{v} ~\mathrm{d} \mathbf{x} + \int_{\Omega} \nabla p \cdot \mathbf{v} ~\mathrm{d} \mathbf{x} = \int_{\Omega} \mathbf{f}\cdot \mathbf{v} ~\mathrm{d} \mathbf{x}
$$

The integration by parts is performed component by component. We impose the essential boundary condition $\mathbf{v}=0$ on $\partial\Omega$, and we denote by

$$
\int_{\Omega}\nabla \mathbf{u}:\nabla \mathbf{v} ~\mathrm{d} \mathbf{x} =\sum_{i=1}^3 \int_{\Omega}\nabla u_i\cdot\nabla v_i ~\mathrm{d} \mathbf{x} =\sum_{i,j=1}^3 \int_{\Omega}\partial_j u_i\partial_j v_i ~\mathrm{d} \mathbf{x} .
$$

We now need to deal with the constraint $\nabla\cdot \mathbf{u}=0$. The theoretical framework for saddle point problems requires that the corresponding bilinear form is the same as the second one appearing in the first part of the variational formulation. To this aim we multiply  $\nabla\cdot \mathbf{u}=0$ by a scalar test function (which will be associated to $p$) and integrate on the full domain, with an integration by parts in order to get the same bilinear form as in the first equation

$$
\int_{\Omega} \nabla\cdot \mathbf{u} ~ q ~\mathrm{d} \mathbf{x}= - \int_{\Omega} \mathbf{u}\cdot \nabla q~\mathrm{d} \mathbf{x}=0,
$$

using that $q=0$ on the boundary as an essential boundary condition.
We finally obtain the mixed variational formulation: 

Find $( \mathbf{u},p)\in H^1_0(\Omega)^3\times H^1_0(\Omega)$ such that

$$
\begin{align}
\left\{
  \begin{array}{llll}
    \int_{\Omega}\nabla \mathbf{u}:\nabla \mathbf{v} ~\mathrm{d} \mathbf{x} &+ \int_{\Omega} \nabla p \cdot \mathbf{v} ~\mathrm{d} \mathbf{x}
    &= \int_{\Omega} \mathbf{f}\cdot \mathbf{v} ~\mathrm{d} \mathbf{x}, & \forall \mathbf{v}\in H_0^1(\Omega)^3 \\
    \int_{\Omega} \mathbf{u}\cdot \nabla q~\mathrm{d} \mathbf{x} & &=0, & \forall p\in H^1_0(\Omega)
  \end{array}
  \right.
\end{align}
$$

#### Examples

- [Stokes problem with mixed (first) formulation on a square domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_stokes_mixed_v1_dir0.py)

#### Second mixed formulation of the Stokes problem

Another possibility to obtained a well posed variational formulation, is to integrate by parts the
$\int_{\Omega} \nabla p \cdot \mathbf{v} ~\mathrm{d} \mathbf{x}$ term in the first formulation:

$$ 
\int_{\Omega} \nabla p \cdot \mathbf{v} ~\mathrm{d} \mathbf{x} = - \int_{\Omega} p \nabla \cdot \mathbf{v} ~\mathrm{d} \mathbf{x} + \int_{\partial\Omega} p \mathbf{v} \cdot \mathbf{n}~\mathrm{d} \sigma= -\int_{\Omega} p ~ \nabla \cdot \mathbf{v} ~\mathrm{d} \mathbf{x} ,
$$

 using here $p=0$ as a natural boundary condition. Note that in the other variational formulation the same boundary condition was essential. In this case, for the second variational formulation, we just multiply $\nabla\cdot \mathbf{u}=0$ by $q$ and integrate. No integration by parts is needed in this case.

$$
\int_{\Omega} \nabla \cdot \mathbf{u} ~ q ~\mathrm{d} \mathbf{x} =0.
$$

This then leads to the following variational formulation:

Find $( \mathbf{u},p)\in H^1(\Omega)^3\times L^2(\Omega)$ such that 

$$
\begin{align}
\left\{
  \begin{array}{llll}
    \int_{\Omega}\nabla \mathbf{u}:\nabla \mathbf{v} ~\mathrm{d} \mathbf{x} &- \int_{\Omega}  p ~ \nabla\cdot \mathbf{v} ~\mathrm{d} \mathbf{x}
    &= \int_{\Omega} \mathbf{f}\cdot \mathbf{v} ~\mathrm{d} \mathbf{x}, &\forall \mathbf{v}\in H^1(\Omega)^3
    \\
    \int_{\Omega}  \nabla\cdot\mathbf{u} ~ q~\mathrm{d} \mathbf{x} & &=0,  &\forall q\in L^2(\Omega)
  \end{array}
  \right.
\end{align}
$$

#### Examples

- [Stokes problem with mixed (second) formulation on a square domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_stokes_mixed_v2_dir0.py)

<a id="linear-elliptic-curl"></a>
### Elliptic-curl Problem

Let $\Omega \subset \mathbb{R}^d$ be an open Liptschitz bounded set, and we look for the solution of the following problem

$$
\begin{align}
  \left\{ 
  \begin{array}{rl}
    \nabla \times \nabla \times \mathbf{u} + \mu \mathbf{u} &= \mathbf{f}, \quad \Omega 
    \\
    \mathbf{u} \times \mathbf{n} &= 0, \quad \partial\Omega
  \end{array} \right.
\end{align}
$$

where $\mathbf{f} \in \mathbf{L}^2(\Omega)$,  $\mu \in L^\infty(\Omega)$ and there exists $\mu_0 > 0$ such that $\mu \geq \mu_0$ almost everywhere.
We take the Hilbert space $V := \mathbf{H}\_0(\mbox{curl}, \Omega)$, in which case the variational formulation writes 

Find $\mathbf{u} \in V$ such that

$$
\begin{align}
  a(\mathbf{u},\mathbf{v}) = l(\mathbf{v}) \quad \forall \mathbf{v} \in V 
\end{align}
where 
\begin{align}
\left\{ 
\begin{array}{rll}
a(\mathbf{u}, \mathbf{v}) &:= \int_{\Omega} \nabla \times \mathbf{u} \cdot \nabla \times \mathbf{v} + \int_{\Omega} \mu \mathbf{u} \cdot \mathbf{v}, & \forall \mathbf{u}, \mathbf{v} \in V  \\
l(\mathbf{v}) &:= \int_{\Omega} \mathbf{v} \cdot \mathbf{f}, & \forall \mathbf{v} \in V  
\end{array} \right.
\end{align}
$$

where 

$$
\begin{align}
\left\{ 
\begin{array}{rll}
a(\mathbf{u}, \mathbf{v}) &:= \int_{\Omega} \nabla \times \mathbf{u} \cdot \nabla \times \mathbf{v} + \int_{\Omega} \mu \mathbf{u} \cdot \mathbf{v}, & \forall \mathbf{u}, \mathbf{v} \in V  \\
l(\mathbf{v}) &:= \int_{\Omega} \mathbf{v} \cdot \mathbf{f}, & \forall \mathbf{v} \in V  
\end{array} \right.
\end{align}
$$

We recall that in $\mathbf{H}\_0(\mbox{curl}, \Omega)$, the bilinear form $a$ is equivalent to the inner product and is therefor continuous and coercive. Hence, our abstract theory applies and there exists a unique solution to the problem.

#### Examples

- [Elliptic-curl problem on a cubic domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/3d_elliptic_curl_dir0.py)

<a id="linear-elliptic-div"></a>
### Elliptic-div Problem

Let $\Omega \subset \mathbb{R}^d$ be an open Liptschitz bounded set, and we look for the solution of the following problem

$$
\begin{align}
  \left\{ 
  \begin{array}{rl}
    - \nabla \nabla \cdot \mathbf{u} + \mu \mathbf{u} &= \mathbf{f}, \quad \Omega 
    \\
    \mathbf{u} \times \mathbf{n} &= 0, \quad \partial\Omega
  \end{array} \right.
\end{align}
$$

where $\mathbf{f} \in \mathbf{L}^2(\Omega)$,  $\mu \in L^\infty(\Omega)$ and there exists $\mu_0 > 0$ such that $\mu \geq \mu_0$ almost everywhere.
We take the Hilbert space $V := \mathbf{H}\_0(\mbox{div}, \Omega)$, in which case the variational formulation writes 


Find $\mathbf{u} \in V$ such that

$$
\begin{align}
  a(\mathbf{u},\mathbf{v}) = l(\mathbf{v}) \quad \forall \mathbf{v} \in V 
\end{align}
$$

where 

$$
\begin{align}
\left\{ 
\begin{array}{rll}
a(\mathbf{u}, \mathbf{v}) &:= \int_{\Omega} \nabla \cdot \mathbf{u} ~ \nabla \cdot \mathbf{v} + \int_{\Omega} \mu \mathbf{u} \cdot \mathbf{v}, & \forall \mathbf{u}, \mathbf{v} \in V  \\
l(\mathbf{v}) &:= \int_{\Omega} \mathbf{v} \cdot \mathbf{f}, & \forall \mathbf{v} \in V  
\end{array} \right.
\end{align}
$$

We recall that in $\mathbf{H}\_0(\mbox{div}, \Omega)$, the bilinear form $a$ is equivalent to the inner product and is therefor continuous and coercive. Hence, our abstract theory applies and there exists a unique solution.

#### Examples

- [Elliptic-div problem on a cubic domain with Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/3d_elliptic_div_dir0.py)

<a id="linear-elasticity"></a>
### Linear Elasticity Problem 

Analysis of deformable structures is essential in engineering, with the equations of linear elasticity being fundamental in this field. In this section, we present the variational formulation of linear elasticity equations using the principle of virtual work.

#### The PDE problem
The governing equations for small elastic deformations of a body \( \Omega \) can be expressed as:
\[
\begin{align}
    -\nabla \cdot \sigma(u) &= f & \text{in } \Omega \\
    \sigma(u) &= \lambda \text{tr}(\epsilon(u))I + 2 \mu \epsilon(u)
\end{align}
\]
where \( \sigma \) is the stress tensor, \( f \) represents the body force per unit volume, \( \lambda \) and \( \mu \) are Lamé's elasticity parameters for the material, \( I \) denotes the identity tensor, and \( \epsilon \) is the symmetric strain tensor. The displacement vector field is denoted by \( u \).

By substituting \( \epsilon(u) \) into \( \sigma \), we obtain:
\[
\sigma(u) = \lambda (\nabla \cdot u)I + \mu(\nabla u + (\nabla u)^T)
\]

#### The Variational Formulation
The variational formulation of the linear elasticity equations involves forming the inner product of the PDE with a vector test function \( v \in \hat{V} \) and integrating over the domain \( \Omega \). This yields:
\[
\int_{\Omega} \sigma : \nabla v \, \mathrm{d} x = \int_{\Omega} f \cdot v \, \mathrm{d} x
\]
Integrating the term \( \nabla \cdot \sigma \cdot v \) by parts, considering boundary conditions, we obtain:
\[
\int_{\Omega} \sigma : \nabla v \, \mathrm{d} x = \int_{\Omega} f \cdot v \, \mathrm{d} x + \int_{\partial \Omega_T} T \cdot v \, \mathrm{d} s
\]
where \( T \) represents the traction vector on the part \( \partial \Omega_T \) of the boundary where it's prescribed.

This leads to the variational formulation: Find \( u \in V \) such that
\[
a(u, v) = L(v) \quad \forall v \in \hat{V}
\]
where
\[
\begin{align}
    a(u, v) &= \int_{\Omega} \sigma(u) : \nabla v \, \mathrm{d} x \\
    L(v) &= \int_{\Omega} f \cdot v \, \mathrm{d} x + \int_{\partial \Omega_T} T \cdot v \, \mathrm{d} s
\end{align}
\]
This formulation can be alternatively expressed as:
\[
a(u, v) = \int_{\Omega} \sigma(u) : \epsilon(v) \, \mathrm{d} x
\]
where \( \epsilon(v) = \frac{1}{2} (\nabla v + (\nabla v)^T) \) is the symmetric strain tensor.

This variational formulation is essential for solving linear elasticity problems numerically using methods like the finite element method (FEM).

- [Linear Elasticity problem on a cubic domain](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/3d_linear_elasticity.py)


<a id="nonlinear-problems"></a>
## Nonlinear Problems

<a id="nonlinear-poisson-2d"></a>
### Nonlinear Poisson in 2D

In this section, we consider the non-linear Poisson problem:

$$
\begin{align}
-\nabla \cdot \left( (1+u^2) \nabla u \right) &= f, \Omega
\\
u &= 0, \partial \Omega
\end{align}
$$

where $\Omega$ denotes the unit square.

For testing, we shall take a function $u$ that fulfills the boundary condition, the compute $f$ as

$$
f(x,y) = -\nabla^2 u - F(u)
$$

The weak formulation is

$$
\int_{\Omega} (1+u^2) \nabla u \cdot \nabla v ~ d\Omega = \int_{\Omega} f v ~d\Omega, \quad \forall v \in \mathcal{V}
$$

For the sack of generality, we shall consider the linear form

$$
G(v;u,w) := \int_{\Omega} (1+w^2) \nabla u \cdot \nabla v ~ d\Omega, \quad \forall u,v,w \in \mathcal{V}
$$

Our problem is then

$$
\begin{align}
\mbox{Find } u \in \mathcal{V}, \mbox{such that}\\
G(v;u,u) = l(v), \quad \forall v \in \mathcal{V}
\end{align}
$$

where

$$
l(v) := \int_{\Omega} f v ~d\Omega, \quad \forall v \in \mathcal{V}
$$

As usual, we'll follow the classical steps,

1. Create a domain
2. Create the Function space
3. Define the different linear/bilinear forms
4. Define the Picard/Newton iteration

#### 1. The topological domain
```python
domain = Square()
B_dirichlet_0 = domain.boundary
```
#### 2. Function space 
```python
V  = ScalarFunctionSpace('V', domain)
```
#### 3.1. Define the Linear form $G$
```python
u  = element_of(V, name='u')
v  = element_of(V, name='v')
w  = element_of(V, name='w')

# Linear form g: V --> R
g = LinearForm(v, integral(domain, (1+w**2)*dot(grad(u), grad(v))))
```
#### 3.2. Define the linear form $L$ 
```python
solution = sin(pi*x)*sin(pi*y)
f = 2*pi**2*(sin(pi*x)**2*sin(pi*y)**2 + 1)*sin(pi*x)*sin(pi*y) - 2*pi**2*sin(pi*x)**3*sin(pi*y)*cos(pi*y)**2 - 2*pi**2*sin(pi*x)*sin(pi*y)**3*cos(pi*x)**2

# Linear form l: V --> R
l = LinearForm(v, integral(domain, f * v))
```
#### 4.1. Picard iteration

$$
\begin{align}
\mbox{Find } u_{n+1} \in \mathcal{V}\_h, \mbox{such that}\\ 
G(v;u_{n+1},u_n) = l(v), \quad \forall v \in \mathcal{V}\_h
\end{align}
$$

The Picard iteration for our problem can be created this way

```python
un  = element_of(V, name='un')

# Bilinear form a: V x V --> R
a = BilinearForm((u, v), g(v, u=u,w=un))

# Dirichlet boundary conditions
bc = [EssentialBC(u, 0, B_dirichlet_0)]

# Variational problem
equation   = find(u, forall=v, lhs=a(u, v), rhs=l(v), bc=bc)

# Error norms
error  = u - solution
l2norm = Norm(error, domain, kind='l2')
```
#### 4.2. Newton iteration
Let's define 

$$
F(v;u) := G(v;u,u) -l(v), \quad \forall v \in \mathcal{V}
$$

Newton method writes

$$
\begin{align}
\mbox{Find } u_{n+1} \in \mathcal{V}\_h, \mbox{such that} \\
F^{\prime}(\delta u,v; u_n) = - F(v;u_n), \quad \forall v \in \mathcal{V}, \\ 
u_{n+1} := u_{n} + \delta u, \quad \delta u \in \mathcal{V}
\end{align}
$$

##### Computing $F^{\prime}$ the derivative of $F$

**SymPDE** allows you to linearize a linear form and get a bilinear form, using the function **linearize**

```python
F = LinearForm(v, g(v,w=u)-l(v))
du  = element_of(V, name='du')

Fprime = linearize(F, u, trials=du)
```

The Newton iteration is then defined as follows,

```python
# Dirichlet boundary conditions
bc = [EssentialBC(du, 0, B_dirichlet_0)]

# Variational problem
equation   = find(du, forall=v, lhs=Fprime(du, v,u=un), rhs=-F(v,u=un), bc=bc)
```

<a id="nonlinear-burgers-1d"></a>
### 1D Burgers equation

We consider the 1d Burgers equation

$$
\partial_t u + u \partial_x u = \nu \frac{\partial ^2u}{\partial x^2}
$$

$u_0(x) := u(x,t)$ denotes the initial condition.
We choose homogeneous neumann boundary conditions in this example, i.e.
$$\partial_n u = 0, \partial \Omega$$ with $\Omega = (0,1)$

#### Time scheme
We shall use a $\theta$-scheme in this case and consider the following problem

$$
\begin{align}
\frac{u^{n+1}-u^n}{\Delta t} + 
\theta~ u^{n+1} \partial_x u^{n+1} + (1-\theta)~ u^n \partial_x u^n = \theta~\nu \frac{\partial ^2u^{n+1}}{\partial x^2} + (1-\theta)~\nu \frac{\partial ^2u^{n}}{\partial x^2}
\end{align}
$$

hence

$$
\begin{align}
u^{n+1} + \Delta t ~ \theta~ u^{n+1} \partial_x u^{n+1} - \Delta t ~ \theta~\nu \frac{\partial ^2u^{n+1}}{\partial x^2} =
u^{n} - \Delta t ~ (1-\theta)~ u^{n} \partial_x u^{n} + \Delta t ~ (1-\theta)~\nu \frac{\partial ^2u^{n}}{\partial x^2}
\end{align}
$$

from now on, we shall denote by $f^n$ the right hand side of the previous equation

$$f^n := u^{n} - \Delta t ~ (1-\theta)~ u^{n} \partial_x u^{n} + \Delta t ~ (1-\theta)~\nu \frac{\partial ^2u^{n}}{\partial x^2}$$

#### Weak formulation

Let $v \in \mathcal{V}$ be a function test, we have by integrating by parts the highest order term:

$$
\begin{align}
\langle v, u^{n+1}  \rangle + \Delta t ~ \theta~ \langle v, u^{n+1} \partial_x u^{n+1}  \rangle + \Delta t ~ \theta~\nu \langle \frac{\partial v}{\partial x}, \frac{\partial u^{n+1}}{\partial x}  \rangle = \langle v, f^n \rangle
\end{align}
$$

The previous weak formulation is still nonlinear with respect to $u^{n+1}$. We shall then follow the same strategy as for the previous chapter on nonlinear Poisson problem.

The strategy is to define the left hand side as a **LinearForm** with respect to $v$, then linearize it around $u^{n+1}$. We therefor can use either Picard or Newton method to treat the nonlinearity.

We consider the following linear form

$$
\begin{align}
G(v;u,w) := \langle v, u  \rangle + \Delta t ~ \theta~ \langle v, w \partial_x u  \rangle + \Delta t ~ \theta~\nu \langle \frac{\partial v}{\partial x}, \frac{\partial u}{\partial x}  \rangle , \quad \forall u,v,w \in \mathcal{V}
\end{align}
$$

Our problem is then

$$
\begin{align}
\mbox{Find } u^{n+1} \in \mathcal{V}, \mbox{such that}\\
G(v;u^{n+1},u^{n+1}) = l(v), \quad \forall v \in \mathcal{V}
\end{align}
$$

where

$$
l(v) := \int_{\Omega} f^n v ~d\Omega, \quad \forall v \in \mathcal{V}
$$

#### SymPDE code

```python
domain = Line()

V = ScalarFunctionSpace('V', domain)
u  = element_of(V, name='u')
v  = element_of(V, name='v')
w  = element_of(V, name='w')
un  = element_of(V, name='un') # time iteration
uk  = element_of(V, name='uk') # nonlinear solver iteration

x = domain.coordinates

nu = Constant('nu')
theta = Constant('theta')
dt = Constant('dt')
```

#### Defining the Linear form $G$

```python
# Linear form g: V --> R
expr = v * u + dt*theta*v*w*dx(u) + dt*theta*nu*dx(v)*dx(u)
g = LinearForm(v, integral(domain, expr))
```

#### Defining the Linear form $l$

```python
# Linear form l: V --> R
expr = v * un - dt*theta*v*un*dx(un) - dt*theta*nu*dx(v)*dx(un)
l = LinearForm(v, integral(domain, expr))
```


#### Picard Method

$$
\begin{align}
\mbox{Find } u_{n+1} \in \mathcal{V}\_h, \mbox{such that}\\
G(v;u_{n+1},u_n) = l(v), \quad \forall v \in \mathcal{V}\_h
\end{align}
$$

##### Picard iteration

```python
# Variational problem
picard = find(u, forall=v, lhs=g(v, u=u,w=uk), rhs=l(v))
```

#### Newton Method

Let's define 

$$
F(v;u) := G(v;u,u) -l(v), \quad \forall v \in \mathcal{V}
$$

Newton method writes

$$
\begin{align}
\mbox{Find } u_{n+1} \in \mathcal{V}\_h, \mbox{such that}\\
F^{\prime}(\delta u,v; u_n) = - F(v;u_n), \quad \forall v \in \mathcal{V} \\
u_{n+1} := u_{n} + \delta u, \quad \delta u \in \mathcal{V}
\end{align}
$$

##### Newton iteration

```python
F = LinearForm(v, g(v,w=u)-l(v))
du  = element_of(V, name='du')

Fprime = linearize(F, u, trials=du)

# Variational problem
newton = find(du, forall=v, lhs=Fprime(du, v,u=uk), rhs=-F(v,u=uk))
```

<a id="nonlinear-navier-stokes-streamfunction-velocity"></a>
### The streamfunction-velocity formulation of the steady-state Navier-Stokes equations for incompressible fluids

The steady-state Navier Stokes problem for an incompressible fluid, with homogeneous Dirichlet boundary conditions (``no slip'' condition), is defined as

$$
\begin{equation}
\left\{
\begin{aligned}
        \left( \mathbf{u} \cdot \nabla \right) \mathbf{u} + \nabla p - 2 \nabla \cdot \left( \frac{1}{R_e} D\left( \mathbf{u} \right) \right) &= \mathbf{f} && \text{in $\Omega$},
  \\
        \nabla \cdot \mathbf{u} &= 0 && \text{in $\Omega$},
  \\
        \mathbf{u} &= 0 && \text{on $\partial\Omega$},
\end{aligned}
\right.
\end{equation}
$$

where $R_e$ is the Reynolds number and $D\left( \mathbf{u} \right) := \frac{1}{2}\left( \nabla \mathbf{u} + \nabla \mathbf{u}^T \right)$ is the strain rate tensor.
The weak formulation reads

$$
\begin{align}
  \text{find $\mathbf{u} \in W$ such that} \quad 
  a(\mathbf{u},\mathbf{v}) + b(\mathbf{u},\mathbf{v};\mathbf{u}) = l(\mathbf{v}) \quad
  \forall \mathbf{v} \in W
\end{align}
$$

where $W \subset \{ \mathbf{v} \in \mathbf{H}_0^1(\Omega): \nabla \cdot \mathbf{v} = 0 \}$ and

$$
\begin{align*}
  a(\mathbf{u}, \mathbf{v}) := \frac{2}{R_e} \int_{\Omega} D\left( \mathbf{u} \right) : D\left(\mathbf{v} \right) ~d\Omega,
  \quad
  b(\mathbf{u}, \mathbf{v}; \mathbf{w}) := \int_{\Omega} \left( \left( \mathbf{w} \cdot \nabla \right) \mathbf{u} \right) \cdot \mathbf{v} ~d\Omega,
  \quad
  l(\mathbf{v}) := \int_{\Omega} \mathbf{f} \cdot \mathbf{v} ~d\Omega.
\end{align*}
$$

When $\Omega$ is a simply connected 2D domain, there exists a unique function $\psi$ such that $\mathbf{u} = \boldsymbol{\nabla} \times \psi:= \left( \partial_y \psi, - \partial_x \psi \right)$;
substituting this expression for $\mathbf{u}$ into \eqref{eq:steady-navier-stokes} leads to the so-called ``streamfunction-velocity formulation'' of the steady-state Navier-Stokes equations for an incompressible fluid.
For more details we refer the reader to \cite{TAGLIABUE2014277}.
The variational formulation is then given by:

$$
\begin{align}
  \text{find $\psi \in V$ such that} \quad 
  a(\psi,\phi) + b(\psi,\phi;\psi) = l(\phi)
  \quad \forall \phi \in V,
\end{align}
$$

where $V \subset H^2(\Omega)$ and  

$$
\begin{align*}
  a(\psi, \phi) := \int_{\Omega} D\left( \boldsymbol{\nabla} \times \psi \right) : D\left(\boldsymbol{\nabla} \times \phi \right) ~d\Omega,  
  \quad
  b(\psi, \phi; \xi) := \int_{\Omega} \left( \left(\boldsymbol{\nabla} \times \xi \cdot \nabla \right) \boldsymbol{\nabla} \times \psi \right) \cdot \boldsymbol{\nabla} \times \phi ~d\Omega,
  \quad
  l(\phi) := \int_{\Omega} \mathbf{f} \cdot \boldsymbol{\nabla} \times \phi ~d\Omega.
\end{align*}
$$

%As in the previous example, one can use a Picard method, by considering $\xi$ as a \texttt{Field} in Eq. \ref{eq:steady-streamfunction-velocity-vf-a}, which will lead to the Python code described in \ref{code:steady-streamfunction-velocity-picard}. %, or consider $a$ as a \texttt{LinearForm} where $\xi = \psi$ are \texttt{Field} objects, and use the \texttt{NewtonIteration} constructor to get a \texttt{BilinearForm}, as described in the Python code \ref{code:steady-streamfunction-velocity-newton}.


#### Linearization of the streamfunction-velocity linear form

In the sequel, we use SymPDE to linearize the streamfunction-velocity formulation, starting from a \texttt{LinearForm}.
As presented in the previous example, we shall write our weak formulation as a linear form where the velocity $\mathbf{u}$ (\textit{i.e.} \texttt{u}) is considered as a \textit{field}, which is implemented in Python code \ref{code:steady-navier-stokes-nl}.
Using the function \texttt{linearize} allows us to construct the bilinear form associated to the linearization of the given linear form, around the field \texttt{psi}, see Python code \ref{code:steady-streamfunction-velocity-lin}.
To get an \texttt{Equation}, the user can add the left-hand-side for the Navier-Stokes equation to the linear form, the appropriate essential boundary conditions, then use the \texttt{NewtonIteration} as presented in the previous example.


<a id="subdomains"></a>
## Subdomains 

<a id="subdomains-linear-poisson"></a>
### The Poisson problem using Nitsche method

In this section, we consider a domain $\Omega$ which a union of multiple subdomain, *i.e.*

$$
\Omega := \bigcup_{i=1}^n \Omega_i, \quad \Omega_i \bigcap \Omega_j = \emptyset, ~ i \neq j  
$$

We shall denote by $\mathcal{I}$ the set of all internal interfaces of $\Omega$.

We shall also need the following operators

- The jump of the function $u$, defined as $[\\![ u ]\\!] := u|\_{\Omega_{i_1}} -  u|\_{\Omega_{i_2}}$ for two adjacent subdomains $\Omega_{i_1}$ and $\Omega_{i_2}$
- The average of the function $u$, defined as $\{u\} := \frac{1}{2} \left( u|\_{\Omega_{i_1}} +  u|\_{\Omega_{i_2}} \right)$ for two adjacent subdomains $\Omega_{i_1}$ and $\Omega_{i_2}$

We consider the Poisson equation

$$
\begin{align}
  \- \nabla^2 u = f \quad \text{in~$\Omega$}, \quad \quad 
  u = 0            \quad \text{on~$\partial \Omega$}, \quad \quad 
   [\\![ u ]\\!] = 0 \quad \text{on~$\mathcal{I}$}, \quad \quad
   [\\![\partial_n u ]\\!] = 0 \quad \text{on~$\mathcal{I}$}.
\end{align}
$$

The variational formulation reads

$$
\begin{align}
  \text{find $u \in V$ such that} \quad a(u,v) = l(v) \quad \forall v \in V,
\end{align}
$$

where 

- $V \subset H^1(\Omega)$, 
- $a(u,v) := \int_{\Omega} \nabla u \cdot \nabla v ~ d\Omega + \int_{\mathcal{I}} \left( \kappa [\\![ u ]\\!] ~ [\\![ v ]\\!] - [\\![ u ]\\!] ~ \partial_n v - \partial_n u ~ [\\![ v ]\\!] \right) ~ d\mathcal{I} $, 
- $l(v) := \int_{\Omega} f v ~ d\Omega ~ d\Gamma$.

#### Examples

- [Poisson on a square with two subdomains and Homogeneous Boundary Conditions](https://github.com/pyccel/sympde/blob/devel-documentation/docs/examples/2d_poisson_subdomains_nitsche_dir0.py)

