# Algebraic and differential operators evalution rules

## Evaluation of grad operator 

Rule ID | symbolic expression | evaluation 
--- | --- | ---
[OG1] | `grad(f+g)` | `grad(f) + grad(g)` 
[OG2] | `grad(\alpha ~ f)` | `\alpha~grad(f)` 
[OG3] | `grad(f ~ g)` | `f~grad(g) + g~grad(f)` 
[OG4] | `grad(\frac{f}{g})` | `-\frac{f}{g^2}~grad(g) + \frac{1}{g}~grad(f)` 
[OG5] | `grad(F+G)` | `grad(F) + grad(G)` 
[OG6] | `grad(\alpha ~ F)` | `\alpha~grad(F)` 
[OG7] | `grad(dot(F, G))` | `convect(F, G) + convect(G, F) + cross(F, curl(G)) - cross(curl(F), G)` 

