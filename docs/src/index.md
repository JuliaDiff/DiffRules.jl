# DiffRules

```@meta
CurrentModule = DiffRules
```

Many differentiation methods rely on the notion of "primitive" differentiation rules that
can be composed via various formulations of the chain rule. Using DiffRules, you can define
new differentiation rules, query whether or not a given rule exists, and symbolically apply
rules to simple Julia expressions.

Note that DiffRules is *not* a fully-fledged symbolic differentiation tool. It is a (very)
simple global database of common derivative definitions, and was developed with the goal of
improving derivative coverage in downstream tools.

Rules for SpecialFunctions, NaNMath and LogExpFunctions live in package extensions, so
DiffRules itself depends on none of them. Load the package you need alongside DiffRules to
get its rules:

```julia
using DiffRules, SpecialFunctions

DiffRules.hasdiffrule(SpecialFunctions.erf, 1)  # true
```

Without `using SpecialFunctions`, that query returns `false` and `diffrules()` does not list
its rules. Packages generating code from `diffrules()` should therefore load the packages
whose rules they want, and skip rules for modules they do not have in scope.

```@docs
DiffRules.@define_diffrule
DiffRules.diffrule
DiffRules.hasdiffrule
DiffRules.diffrules
```
