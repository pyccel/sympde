# Algebraic and differential operators evalution rules

## Evaluation of grad operator 

Rule ID | symbolic expression | evaluation 
--- | --- | ---
`[OG1]` | $\texttt{grad}(f+g)$ | $\texttt{grad}(f) + \texttt{grad}(g)$ 
`[OG2]` | $\texttt{grad}(\alpha ~ f)$ | $\alpha~\texttt{grad}(f)$ 
`[OG3]` | $\texttt{grad}(f ~ g)$ | $f~\texttt{grad}(g) + g~\texttt{grad}(f)$ 
`[OG4]` | $\texttt{grad}(\frac{f}{g})$ | $-\frac{f}{g^2}~\texttt{grad}(g) + \frac{1}{g}~\texttt{grad}(f)$ 
`[OG5]` | $\texttt{grad}(F+G)$ | $\texttt{grad}(F) + \texttt{grad}(G)$ 
`[OG6]` | $\texttt{grad}(\alpha ~ F)$ | $\alpha~\texttt{grad}(F)$ 
`[OG7]` | $\texttt{grad}(\texttt{dot}(F, G))$ | $\texttt{convect}(F, G) + \texttt{convect}(G, F) + \texttt{cross}(F, \symcurl(G)) - \texttt{cross}(\symcurl(F), G)$ 
