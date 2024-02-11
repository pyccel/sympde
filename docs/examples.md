# Examples

## Linear Problems

### Mixed FEM Poisson

Let $\Omega \subset \mathbb{R}^3$ and consider the Poisson problem

$$
\begin{align}
  \left\{ 
  \begin{array}{clr}
    -\Delta p & =f & ,~\Omega    \\
    p         & =0 & ,~\partial \Omega
  \end{array} \right.
  \label{eq:abs_poisson_mixed}
\end{align}
$$

Using that $\Delta p = \nabla\cdot\nabla p$, we set $ \mathbf{u}=\nabla p$, then the Poisson equation \ref{eq:abs_poisson_mixed} can be written equivalently

$$ \mathbf{u}=-\nabla p, ~~~ \nabla\cdot \mathbf{u}= f.$$

#### First mixed formulation of the Poisson problem

Instead of having one unknown, we now have two, along with the above two equations.
In order to get a mixed variational formulation, we first take the dot product of the first one by $ \mathbf{v}$ and integrate by parts

$$
\begin{align}
\int_{\Omega} \mathbf{u}\cdot \mathbf{v}~\mathrm{d} \mathbf{x} -\int_{\Omega} p\,\nabla\cdot \mathbf{v}~\mathrm{d} \mathbf{x} + \int_{\partial\Omega} p \,  \mathbf{v}\cdot \mathbf{n}~\mathrm{d} \sigma = \int_{\Omega} \mathbf{u}\cdot \mathbf{v}~\mathrm{d} \mathbf{x} -\int_{\Omega} p\,\nabla\cdot \mathbf{v}~\mathrm{d} \mathbf{x}=0,
\end{align}
$$

using $p=0$ as a natural boundary condition. Then multiplying the second equation by $q$ and integrating yields

$$
\begin{align}
\int_{\Omega} \nabla\cdot\mathbf{u} \, q ~\mathrm{d} \mathbf{x} = \int_{\Omega} f q ~\mathrm{d} \mathbf{x}.
\end{align}
$$

No integration by parts is necessary here. And we thus get the following mixed variational formulation:

Find $(\mathbf{u},p) \in H(\mathrm{div},\Omega)\times L^2(\Omega)$ such that

$$
\begin{align}
\left\{ 
\begin{array}{llll}
  \int_{\Omega} \mathbf{u}\cdot \mathbf{v}~\mathrm{d} \mathbf{x} &- \int_{\Omega} p\,\nabla\cdot \mathbf{v}~\mathrm{d} \mathbf{x} &=0, & \forall \mathbf{v}\in H(\mathrm{div},\Omega) \\
  - \int_{\Omega} \nabla\cdot\mathbf{u} \, q ~\mathrm{d} \mathbf{x} &  &= - \int_{\Omega} f q ~\mathrm{d} \mathbf{x}, & \forall q\in L^2(\Omega)
\end{array} \right.
\end{align}
$$

todo

```python
# ... abstract model
domain = Square()

V1 = VectorFunctionSpace('V1', domain, kind='Hdiv')
V2 = ScalarFunctionSpace('V2', domain, kind='L2')
X  = ProductSpace(V1, V2)

x,y = domain.coordinates

# rhs
f0 = -2*x*(1-x) -2*y*(1-y)
# analytical solution
u  = x*(1-x)*y*(1-y)

F = element_of(V2, name='F')

p,q = [element_of(V1, name=i) for i in ['p', 'q']]
u,v = [element_of(V2, name=i) for i in ['u', 'v']]

int_0 = lambda expr: integral(domain , expr)

a  = BilinearForm(((p,u),(q,v)), int_0(dot(p,q) + div(q)*u + div(p)*v) )
l  = LinearForm((q,v), int_0(f0*v))

# ...
error = F-u
l2norm_F = Norm(error, domain, kind='l2')

# ...
equation = find([p,u], forall=[q,v], lhs=a((p,u),(q,v)), rhs=l(q,v))
```

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
\int_{\Omega}\nabla \cdot\mathbf{u}  \, q ~\mathrm{d} \mathbf{x} = -\int_{\Omega}  \mathbf{u} \cdot \nabla q  ~\mathrm{d} \mathbf{x} = \int_{\Omega} f q ~\mathrm{d} \mathbf{x},
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

```python
# TODO
```


## Nonlinear Problems

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
