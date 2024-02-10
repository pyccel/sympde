# Algebraic and differential operators evalution rules

## Evaluation of the $\mathrm{grad}$ operator 

Rule ID | symbolic expression | evaluation 
--- | --- | ---
`[OG1]` | $\mathrm{grad}(f+g)$ | $\mathrm{grad}(f) + \mathrm{grad}(g)$ 
`[OG2]` | $\mathrm{grad}(\alpha ~ f)$ | $\alpha~\mathrm{grad}(f)$ 
`[OG3]` | $\mathrm{grad}(f ~ g)$ | $f~\mathrm{grad}(g) + g~\mathrm{grad}(f)$ 
`[OG4]` | $\mathrm{grad}(\frac{f}{g})$ | $-\frac{f}{g^2}\mathrm{grad}(g) + \frac{1}{g}\mathrm{grad}(f)$ 
`[OG5]` | $\mathrm{grad}(F+G)$ | $\mathrm{grad}(F) + \mathrm{grad}(G)$ 
`[OG6]` | $\mathrm{grad}(\alpha ~ F)$ | $\alpha~\mathrm{grad}(F)$ 
`[OG7]` | $\mathrm{grad}(\mathrm{dot}(F, G))$ | $\mathrm{convect}(F, G) + \mathrm{convect}(G, F) + \mathrm{cross}(F, \mathrm{curl}(G)) - \mathrm{cross}(\mathrm{curl}(F), G)$ 

## Evaluation of $\mathrm{curl}$ and $\mathrm{rot}$ operators in 2D 

Rule ID | symbolic expression | evaluation 
--- | --- | ---
`[OC1a]` | $\mathrm{curl}(f+g)$ | $\mathrm{curl}(f) + \mathrm{curl}(g)$ 
`[OC2a]` | $\mathrm{curl}(\alpha ~ f)$ | $\alpha~\mathrm{curl}(f)$ 
`[OC1b]` | $\mathrm{rot}(F+G)$ | $\mathrm{rot}(F) + \mathrm{rot}(G)$ 
`[OC2b]` | $\mathrm{rot}(\alpha ~ F)$ | $\alpha~\mathrm{rot}(F)$ 


## Evaluation of $\mathrm{curl}$ operator in 3D

Rule ID | symbolic expression | evaluation 
--- | --- | ---
`[OC1]` | $\mathrm{curl}(F+G)$ | $\mathrm{curl}(F) + \mathrm{curl}(G)$ 
`[OC2]` | $\mathrm{curl}(\alpha ~ F)$ | $\alpha~\mathrm{curl}(F)$ 
`[OC3]` | $\mathrm{curl}(f  ~ F)$ | $f~\mathrm{curl}(F) + \mathrm{cross}(\mathrm{grad}(f), F)$ 
`[OC4]` | $\mathrm{curl}(\mathrm{cross}(F, G))$ | $\mathrm{div}(G)~F - \mathrm{div}(F)~G - \mathrm{convect}(F, G) +\mathrm{convect}(G, F)$ 

## Evaluation of $\mathrm{div}$ operator 

Rule ID | symbolic expression | evaluation 
--- | --- | ---
`[OD1]` | $\mathrm{div}(F+G)$ | $\mathrm{div}(F) + \mathrm{div}(G)$ 
`[OD2]` | $\mathrm{div}(\alpha ~ F)$ | $\alpha~\mathrm{div}(F)$ 
`[OD3]` | $\mathrm{div}(f ~ G)$ | $f~\mathrm{div}(G) + \mathrm{dot}(G, \mathrm{grad}(f))$ 
`[OD4]` | $\mathrm{div}(\mathrm{cross}(F, G))$ | $-\mathrm{dot}(F, \mathrm{curl}(G)) + \mathrm{dot}(G, \mathrm{curl}(F))$ 

## Evaluation of $\mathrm{laplace}$ operator 

Rule ID | symbolic expression | evaluation 
--- | --- | ---
`[OL1]` | $\mathrm{laplace}(f+g)$ | $\mathrm{laplace}(f) + \mathrm{laplace}(g)$ 
`[OL2]` | $\mathrm{laplace}(\alpha ~ f)$ | $\alpha~\mathrm{laplace}(f)$ 
`[OL3]` | $\mathrm{laplace}(f~g)$ | $f ~ \mathrm{laplace}(g) + g ~ \mathrm{laplace}(f) + 2 \mathrm{dot}(\mathrm{grad}(f), \mathrm{grad}(g))$ 

## Evaluation of specific combination of operators 

Rule ID | symbolic expression | evaluation 
--- | --- | ---
`[OS1]` | $\mathrm{curl}(\mathrm{grad}(f))$ | $0$ 
`[OS2]` | $\mathrm{div}(\mathrm{curl}(F))$ | $0$ 
`[OS3]` | $\mathrm{div}(\mathrm{cross}(\mathrm{grad}(F), \mathrm{grad}(G)))$ | $0$ 
`[OS4]` | $\mathrm{curl}(\mathrm{curl}(F))$ | $\mathrm{grad}(\mathrm{div}(F)) - \mathrm{laplace}(F)$ 
`[OS5]` | $\mathrm{curl}(f~\mathrm{grad}(g))$ | $\mathrm{cross}(\mathrm{grad}(f), \mathrm{grad}(g))$ 

