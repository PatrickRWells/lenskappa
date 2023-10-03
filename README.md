This is a set of tools that are being used in the TDCOSMO collaboration to
study the environments around strong gravitational lenses.

Most of this is not stand-alone code. Weighted number counts are analysis templates
that are designed to be used by [cosmap](https://github.com/PatrickRWells/cosmap).

The code in the "inference" module is used to convert outputs of the `cosmap` analysis
and compute a histogram on $\kappa_{ext}$. It also contains routines for doing
hierarchical modeling of populations of strong lens lines of sight.