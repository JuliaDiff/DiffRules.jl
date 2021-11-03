# DiffRules

[![Github CI Status](https://github.com/JuliaDiff/DiffRules.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/JuliaDiff/DiffRules.jl/actions/workflows/ci.yml)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiff/DiffRules.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaDiff/DiffRules.jl?branch=master)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliadiff.org/DiffRules.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliadiff.org/DiffRules.jl/latest)

Many differentiation methods rely on the notion of "primitive" differentiation rules that
can be composed via various formulations of the chain rule. Using DiffRules, you can define
new differentiation rules, query whether or not a given rule exists, and symbolically apply
rules to simple Julia expressions.

Note that DiffRules is *not* a fully-fledged symbolic differentiation tool. It is a (very)
simple global database of common derivative definitions, and was developed with the goal of
improving derivative coverage in downstream tools.
