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

```@docs
DiffRules.@define_diffrule
DiffRules.diffrule
DiffRules.hasdiffrule
DiffRules.diffrules
```
