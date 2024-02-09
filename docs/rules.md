# Algebraic and differential operators evalution rules

\begin{table}[h!]
  \centering
  {\small
  \begin{tabular}{lll}
    \hline
    Rule ID & symbolic expression & evaluation \\
    \hline
    \texttt{[OG1]} & $\symgrad(\symf+\symg)$ & $\symgrad(\symf) + \symgrad(\symg)$ \\
    \texttt{[OG2]} & $\symgrad(\alpha ~ \symf)$ & $\alpha~\symgrad(\symf)$ \\
    \texttt{[OG3]} & $\symgrad(\symf ~ \symg)$ & $\symf~\symgrad(\symg) + \symg~\symgrad(\symf)$ \\
    \texttt{[OG4]} & $\symgrad(\frac{\symf}{\symg})$ & $-\frac{\symf}{\symg^2}~\symgrad(\symg) + \frac{1}{\symg}~\symgrad(\symf)$ \\
    \texttt{[OG5]} & $\symgrad(\symF+\symG)$ & $\symgrad(\symF) + \symgrad(\symG)$ \\
    \texttt{[OG6]} & $\symgrad(\alpha ~ \symF)$ & $\alpha~\symgrad(\symF)$ \\
    \texttt{[OG7]} & $\symgrad(\symdot(\symF, \symG))$ & $\symconvect(\symF, \symG) + \symconvect(\symG, \symF)$ \\
                   &                   & $+ \symcross(\symF, \symcurl(\symG)) - \symcross(\symcurl(\symF), \symG)$ \\

    \hline
  \end{tabular}
  }
  \caption{Evaluation of $\symgrad$ operator}
  \label{tab:sympde-op-grad}
\end{table}

\begin{table}[h!]
  \centering
  \begin{tabular}{lll}
    \hline
    Rule ID & symbolic expression & evaluation \\
    \hline
    \texttt{[OC1a]} & $\symcurl(\symf+\symg)$ & $\symcurl(\symf) + \symcurl(\symg)$ \\
    \texttt{[OC2a]} & $\symcurl(\alpha ~ \symf)$ & $\alpha~\symcurl(\symf)$ \\
    \texttt{[OC1b]} & $\symrot(\symF+\symG)$ & $\symrot(\symF) + \symrot(\symG)$ \\
    \texttt{[OC2b]} & $\symrot(\alpha ~ \symF)$ & $\alpha~\symrot(\symF)$ \\
    \hline
  \end{tabular}
  \caption{Evaluation of $\symcurl$ and $\symrot$ operators in 2D}
  \label{tab:sympde-op-curl-2d}
\end{table}

\begin{table}[h!]
  \centering
  \begin{tabular}{lll}
    \hline
    Rule ID & symbolic expression & evaluation \\
    \hline
    \texttt{[OC1]} & $\symcurl(\symF+\symG)$ & $\symcurl(\symF) + \symcurl(\symG)$ \\
    \texttt{[OC2]} & $\symcurl(\alpha ~ \symF)$ & $\alpha~\symcurl(\symF)$ \\
    \texttt{[OC3]} & $\symcurl(\symf  ~ \symF)$ & $\symf~\symcurl(\symF) + \symcross(\symgrad(f), \symF)$ \\
    \texttt{[OC4]} & $\symcurl(\symcross(\symF, \symG))$ & $\symdiv(\symG)~\symF - \symdiv(\symF)~\symG$ \\
                   &                     & $- \symconvect(\symF, \symG) +\symconvect(\symG, \symF)$ \\

    \hline
  \end{tabular}
  \caption{Evaluation of $\symcurl$ operator in 3D}
  \label{tab:sympde-op-curl-3d}
\end{table}

\begin{table}[h!]
  \centering
  \begin{tabular}{lll}
    \hline
    Rule ID & symbolic expression & evaluation \\
    \hline
    \texttt{[OD1]} & $\symdiv(\symF+\symG)$ & $\symdiv(\symF) + \symdiv(\symG)$ \\
    \texttt{[OD2]} & $\symdiv(\alpha ~ \symF)$ & $\alpha~\symdiv(\symF)$ \\
    \texttt{[OD3]} & $\symdiv(\symf ~ \symG)$ & $\symf~\symdiv(\symG) + \symdot(\symG, \symgrad(\symf))$ \\
    \texttt{[OD4]} & $\symdiv(\symcross(\symF, \symG))$ & $-\symdot(\symF, \symcurl(\symG)) + \symdot(\symG, \symcurl(\symF))$ \\
    \hline
  \end{tabular}
  \caption{Evaluation of $\symdiv$ operator}
  \label{tab:sympde-op-div}
\end{table}

\begin{table}[h!]
  \centering
  {\small
  \begin{tabular}{lll}
    \hline
    Rule ID & symbolic expression & evaluation \\
    \hline
    \texttt{[OL1]} & $\symlaplace(\symf+\symg)$ & $\symlaplace(\symf) + \symlaplace(\symg)$ \\
    \texttt{[OL2]} & $\symlaplace(\alpha ~ \symf)$ & $\alpha~\symlaplace(\symf)$ \\
    \texttt{[OL3]} & $\symlaplace(\symf~\symg)$ & $\symf ~ \symlaplace(\symg) + \symg ~ \symlaplace(\symf)$ \\
                   &           & $+ 2 \symdot(\symgrad(\symf), \symgrad(\symg))$ \\
    \hline
  \end{tabular}
  }
  \caption{Evaluation of $\symlaplace$ operator}
  \label{tab:sympde-op-laplace}
\end{table}

\begin{table}[h!]
  \centering
  {\small
  \begin{tabular}{lll}
    \hline
    Rule ID & symbolic expression & evaluation \\
    \hline
    \texttt{[OS1]} & $\symcurl(\symgrad(\symf))$ & $0$ \\
    \texttt{[OS2]} & $\symdiv(\symcurl(\symF))$ & $0$ \\
    \texttt{[OS3]} & $\symdiv(\symcross(\symgrad(\symF), \symgrad(\symG)))$ & $0$ \\
    \texttt{[OS4]} & $\symcurl(\symcurl(\symF))$ & $\symgrad(\symdiv(\symF)) - \symlaplace(\symF)$ \\
    \texttt{[OS5]} & $\symcurl(\symf~\symgrad(\symg))$ & $\symcross(\symgrad(\symf), \symgrad(\symg))$ \\
    \hline
  \end{tabular}
  }
  \caption{Evaluation of specific combination of operators}
  \label{tab:sympde-op-prop}
\end{table}

