

\begin{table}[h!]
  \centering
  {\scriptsize
    \begin{tabular}{ll}
      \hline
      sympde notation & Mathematical notion \\
      \hline
      \hline
      \texttt{Domain('Omega', dim=n)}  & an abstract domain $\Omega \subset \mathbb{R}^n$\\
      \texttt{Domain('Omega', interiors, boundaries, connectivity)}  & an abstract domain $\Omega \subset \mathbb{R}^n$,\\
                                                                     & represented by its interiors, boundaries and connectivity\\
      \texttt{Interval('I', coordinate=x)}  & an abstract (interior) interval $I \subset \mathbb{R}$, with the coordinate $x$\\
      \texttt{Union(A, B, C)}  & union of domains $A \cup B \cup C$\\
      \texttt{ProductDomain(A, B)}  & product of two domains $A \times B$\\
      \texttt{Indicator(A)}  & Indicator/characteristic function on $A$, \textit{i.e.} $\mathbf{1}_A$\\
      \texttt{Line('Omega')}    & an abstract line segment, $\Omega \subset \mathbb{R}$, \\
                                & having one interval as interior domain\\
      \texttt{Square('Omega')}  & an abstract square, $\Omega \subset \mathbb{R}^2$, \\
                                & having one interior domain as product of two intervals \\
      \texttt{Cube('Omega')}    & an abstract cube, $\Omega \subset \mathbb{R}^3$, \\
                                & having one interior domain as product of three intervals \\
      \texttt{NormalVector('n')} & Normal vector $\mathbf{n}$ \\
      \texttt{TangentVector('t')} & Tangent vector $\mathbf{t}$ \\
      \texttt{e = ElementDomain(Omega)} & an element of an abstract domain $e \in \mathcal{T}(\Omega)$ or $e = d\Omega$,  \\
                                   & where $\mathcal{T}(\Omega)$ is a given tessellation of $\Omega$ \\
      \texttt{DomainArea(Omega)} & Area of an abstract domain $\Omega$ \\
      \texttt{ElementArea(Omega)} & Area of an abstract element of a domain $\Omega$ \\
      \texttt{Area(A)} & Area of an expression of topological domain notions \\
      \hline
    \end{tabular}
    \caption{Table of related notions to a topological domain and their mathematical meaning}
    \label{tab:sympde-domain}
  }
\end{table}



\begin{table}[h!]
  \centering
  {\scriptsize
    \begin{tabular}{ll}
      \hline
      sympde notation & Mathematical notion \\
      \hline
      \hline
      \texttt{Mapping('F', n)}          & a mapping functor $\mathbf{F}_{\Omega} := \left( \mathbf{F}, \Omega \right): \Omega \rightarrow \mathbb{R}^n  $\\
      \texttt{F[i]}            & $i^{th}$ physical coordinate from $\{x,y,z\}$\\
      \texttt{F.jacobian}      & the jacobian matrix $\mathcal{D}_\mathbf{F}$\\
      \texttt{F.det\_jacobian} & determinant of the jacobian matrix, $J_\mathbf{F}:= \mathrm{det} ~\mathcal{D}_\mathbf{F}$\\
      \texttt{F.covariant}     & the covariant matrix $\left( \mathcal{D}_\mathbf{F} \right)^{-T}$\\
      \texttt{F.contravariant} & the contravariant matrix $\frac{1}{J_\mathbf{F}} \mathcal{D}_\mathbf{F}$ \\
      \texttt{F.hessian}       & the hessian matrix $\mathcal{H}_\mathbf{F}$\\
      \hline
      \hline
      \texttt{Covariant(F, v)} & action of the covariant matrix of $\mathbf{F}$ on $\mathbf{v}$, \textit{i.e.} $\left( \mathcal{D}_\mathbf{F} \right)^{-T} \mathbf{v}$\\
      \texttt{Contravariant(F, v)} & action of the contravariant matrix of $\mathbf{F}$ on $\mathbf{v}$, \textit{i.e.} $\frac{1}{J_\mathbf{F}} \mathcal{D}_\mathbf{F} \mathbf{v}$\\
      \hline
    \end{tabular}
    \caption{Table of related notions to a mapping their mathematical meaning}
    \label{tab:sympde-mapping}
  }
\end{table}



\begin{table}[h!]
  \centering
  {\scriptsize
    \begin{tabular}{ll}
      \hline
      sympde notation & Mathematical notion \\
      \hline
      \hline
      \texttt{ScalarFunctionSpace('V', Omega)}     & scalar function space $\mathcal{V}$\\
      \texttt{VectorFunctionSpace('W', Omega)}     & vector function space $\mathbf{\mathcal{W}}$\\
      \texttt{ProductSpace(V, W)} or \texttt{V*W}  & product of spaces, \textit{i.e.} $\mathcal{V} \times \mathcal{W}$\\
      \hline
      \hline
      \texttt{element\_of\_space(V, 'v')}          & scalar function $v \in \mathcal{V}$ \\
      \texttt{element\_of\_space(W, 'w')}          & vector function $\mathbf{w} \in \mathcal{W}$  \\
      \texttt{element\_of\_space(V*W, ['v', 'w'])} & $\left(v,\mathbf{w}\right) \in \mathcal{V} \times \mathcal{W}$ \\
      \hline
    \end{tabular}
    \caption{Table of function spaces and how to create their elements in SymPDE}
    \label{tab:sympde-concepts}
  }
\end{table}


\begin{table}[h!]
  \centering
  \begin{tabular}{ll}
    \hline
    sympde notation & Mathematical notion \\
    \hline
    \texttt{H1SpaceType}    & $H^1$    \\
    \texttt{HcurlSpaceType} & ${H}{\mbox{curl}}$ \\
    \texttt{HdivSpaceType}  & ${H}{\mbox{div}}$  \\
    \texttt{L2SpaceType}    & $L^2$    \\
    \hline
  \end{tabular}
  \caption{Function space types}
  \label{tab:sympde-space-type}
\end{table}


\begin{table}[h!]
  \centering
  \begin{tabular}{ll}
    \hline
    sympde notation & Mathematical notion \\
    \hline
    \texttt{ScalarFunctionSpace('V0', Omega, kind='H1')}          & scalar function space $\mathcal{V}_0 \subseteq H^1(\Omega)$    \\
    \texttt{VectorFunctionSpace('V1', Omega, kind='Hcurl')} & vector function space $\mathcal{V}_1 \subseteq H(\mbox{curl}, \Omega)$ \\
    \texttt{VectorFunctionSpace('V2', Omega, kind='Hdiv')}  & vector function space $\mathcal{V}_2 \subseteq H(\mbox{div}, \Omega)$  \\
    \texttt{ScalarFunctionSpace('V3', Omega, kind='L2')}          & scalar function space $\mathcal{V}_3 \subseteq L^2(\Omega)$    \\
    \hline
  \end{tabular}
  \caption{Declaration of typed function space in 3D}
  \label{tab:sympde-space-type-declaration}
\end{table}



\begin{table}[h!]
  \centering
    \begin{tabular}{ll}
      \hline
      sympde notation & Mathematical notion \\
      \hline
      \texttt{P\_V := Projector(V)}   & Projector onto the function space $\mathcal{V}$  \\
      \texttt{P\_V(expr)}   & Projection of the expression \texttt{expr} onto the function space $\mathcal{V}$  \\
      \texttt{Pi\_V := Projector(V, 'commuting')}   & Commuting projector onto the typed function space $\mathcal{V}$  \\
      \texttt{Pi\_V(expr)}   & Commuting projection of the expression \texttt{expr} onto the typed function space $\mathcal{V}$  \\
      \hline
    \end{tabular}
    \caption{Table of Projection on function space and their related mathematical notion.}
    \label{tab:projection}
\end{table}


\begin{table}[h!]
  \centering
  \begin{minipage}[t]{0.4\textwidth}
  {\scriptsize
    \begin{tabular}{ll}
      \hline
      Atomic variable & Mathematical meaning \\
      \hline
      \texttt{ScalarFunction(V, 'u')} &  $u \in \mathcal{V}$\\
      \texttt{VectorFunction(W, 'w')} &  $\mathbf{w} \in \mathcal{W}$\\
            \texttt{Constant('mu', real=True)} & constant number $\mu \in \mathbb{R}$ \\
      Python \texttt{int}, \texttt{float} & integer/floating-point numbers \\
      \hline
    \end{tabular}
    \caption{Table of atomic variables.}
    \label{tab:sympde-atomic}
  }
  \end{minipage}
\end{table}




\begin{table}[h!]
  \centering
  \begin{minipage}[t]{0.4\textwidth}
  {\scriptsize
    \begin{tabular}{ll}
      \hline
      Mathematical notion & SymPDE notation \\
      \hline
      $\mathbf{u} \cdot \mathbf{v}$ & \texttt{dot(u, v)} \\
      $\mathbf{u} : \mathbf{v}$ & \texttt{inner(u, v)} \\
      $\mathbf{u} \times \mathbf{v}$ & \texttt{cross(u, v)} \\
      $\mathbf{u} \otimes \mathbf{v}$ & \texttt{outer(u, v)} \\
      \hline
    \end{tabular}
    \caption{Table of algebraic operators.}
    \label{tab:sympde-algebra}
  }
  \end{minipage}
  \hfill
  \begin{minipage}[t]{0.5\textwidth}
  {\scriptsize
    \begin{tabular}{ll}
      \hline
      Mathematical notion & SymPDE notation  \\
      \hline
      $\partial_x u$,~ $\partial_y u$~ or~ $\partial_z u$  & \texttt{dx(u)},~ \texttt{dy(u)}~ or~ \texttt{dz(u)} \\
      $\nabla u$~ or~ $\nabla \mathbf{u}$   & \texttt{grad(u)} \\
      $\nabla \cdot \mathbf{u}$   & \texttt{div(u)} \\
      $\nabla \times \mathbf{u}$   & \texttt{curl(u)} \\
      2D rotational $\nabla \times u$   & \texttt{rot(u)} \\
      $\left( \mathbf{a} \cdot \nabla \right) \mathbf{u}$   & \texttt{convect(a, u)} \\
      trace $\gamma(u)$   & \texttt{trace(u)} \\
      normal derivative $\partial_{\mathbf{n}} u$   & \texttt{Dn(u)} \\
      Strain tensor $\frac{1}{2}\left(\nabla \mathbf{u} + \nabla \mathbf{u}^T \right) $   & \texttt{D(u)} \\
      Laplace $\Delta u$   & \texttt{laplace(u)} \\
      Hessian $H(u)$   & \texttt{hessian(u)} \\
      Poisson Bracket $[u,v]$   & \texttt{bracket(u,v)} \\
      \hline
    \end{tabular}
    \caption{Table of differential operators (boldface font is used for vector functions/expressions).}
    \label{tab:sympde-diffops}
  }
  \end{minipage}
\end{table}

\begin{table}
  \centering
  \begin{minipage}[t]{0.5\textwidth}
  {\small
    \begin{tabular}{ll}
      \hline
      Mathematical notion & SymPDE notation  \\
      \hline
      jump of $u$, \textit{i.e.} $[u]$   & \texttt{jump(u)} \\
      average of $u$, \textit{i.e.} $\langle u \rangle$   & \texttt{avg(u)} \\
      convolution of $u$ with a kernel $K$, \textit{i.e.} $K * u$   & \texttt{conv(K,u)} \\
      \hline
    \end{tabular}
    \caption{Table of other operators.}
    \label{tab:sympde-otherops}
  }
  \end{minipage}
\end{table}



\begin{table}[h!]
  \centering
  {\small
    \begin{tabular}{ll}
      \hline
      sympde notation & Mathematical notion \\
      \hline
      \hline
      \texttt{DomainIntegral(f)}   & integral over a domain, \textit{i.e.} $(f, \Omega) \mapsto \int_{\Omega} f ~d\Omega$\\
      \texttt{BoundaryIntegral(f)} & integral over a boundary, \textit{i.e.} $(f, \Gamma) \mapsto \int_{\Gamma} f ~d\partial\Omega$\\
      \texttt{PathIntegral(F)}     & integral over an oriented path, \textit{i.e.} $(\mathbf{F}, C) \mapsto \int_{C} \mathbf{F} \cdot d\mathbf{s}$\\
      \texttt{SurfaceIntegral(F)}  & integral over an oriented surface, \textit{i.e.} $(\mathbf{F}, S) \mapsto \int_{S} \mathbf{F} \cdot d\mathbf{S}$\\
      \hline
    \end{tabular}
    \caption{Integral operators in SymPDE and their mathematical meaning.}
    \label{tab:integrals}
  }
\end{table}



