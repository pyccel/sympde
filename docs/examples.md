# Examples

## Nonlinear Problems

### Nonlinear Poisson in 2D

In this section, we consider the non-linear Poisson problem:

$$
-\nabla \cdot \left( (1+u^2) \nabla u \right) = f, \Omega
\\
u = 0, \partial \Omega
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
\mbox{Find } u \in \mathcal{V}, \mbox{such that}\\
G(v;u,u) = l(v), \quad \forall v \in \mathcal{V}
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
#### 4. Picard iteration

$$
\mbox{Find } u_{n+1} \in \mathcal{V}\_h, \mbox{such that} 

G(v;u_{n+1},u_n) = l(v), \quad \forall v \in \mathcal{V}\_h
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
#### 5. Newton iteration
Let's define 

$$
F(v;u) := G(v;u,u) -l(v), \quad \forall v \in \mathcal{V}
$$

Newton method writes

$$
\mbox{Find } u_{n+1} \in \mathcal{V}\_h, \mbox{such that}

F^{\prime}(\delta u,v; u_n) = - F(v;u_n), \quad \forall v \in \mathcal{V}, \\ 
u_{n+1} := u_{n} + \delta u, \quad \delta u \in \mathcal{V}
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
