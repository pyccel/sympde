# Algebraic and differential operators evalution rules

## Evaluation of grad operator 

Rule ID | symbolic expression | evaluation 
--- | --- | ---
`[OG1]` | $\mathrm{grad}(f+g)$ | $\mathrm{grad}(f) + \mathrm{grad}(g)$ 
`[OG2]` | $\mathrm{grad}(\alpha ~ f)$ | $\alpha~\mathrm{grad}(f)$ 
`[OG3]` | $\mathrm{grad}(f ~ g)$ | $f~\mathrm{grad}(g) + g~\mathrm{grad}(f)$ 
`[OG4]` | $\mathrm{grad}(\frac{f}{g})$ | $-\frac{f}{g^2}~\mathrm{grad}(g) + \frac{1}{g}~\mathrm{grad}(f)$ 
`[OG5]` | $\mathrm{grad}(F+G)$ | $\mathrm{grad}(F) + \mathrm{grad}(G)$ 
`[OG6]` | $\mathrm{grad}(\alpha ~ F)$ | $\alpha~\mathrm{grad}(F)$ 
`[OG7]` | $\mathrm{grad}(\mathrm{dot}(F, G))$ | $\mathrm{convect}(F, G) + \mathrm{convect}(G, F) + \mathrm{cross}(F, \mathrm{curl}(G)) - \mathrm{cross}(\mathrm{curl}(F), G)$ 
