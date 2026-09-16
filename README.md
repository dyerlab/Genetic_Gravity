# Genetic_Gravity

The latest build of this entire manuscript is found [here](https://dyerlab.github.io/Genetic_Gravity/).

---

## Summary

> Population Graphs summarize genetic connectivity among populations as an undirected network, but the direction of gene flow along any connection — who is a source, who is a sink — is not something the topology can show on its own. This paper introduces *genetic gravity*, a way to recover that missing direction without adding a new demographic model, an equilibrium assumption, spatial coordinates, or any tunable parameter. For each population, a local kernel (scaled to its own neighborhood of genetic distances) turns every symmetric edge into a pair of directional weights; the difference between them, Δ, is an asymmetry index that can be read directly off a single Population Graph. Using individual-based simulations under symmetric and asymmetric migration, we show the index is specific (it stays near zero under symmetric migration no matter how much genetic differentiation has built up), most reliable as a graph-wide average rather than edge by edge, has characterizable detection limits governed by the migration ratio and effective population size, and avoids the drift-driven false positives that affect differentiation-based tools such as divMigrate. The method, along with a hierarchy of significance tests and permutation nulls for applying it to real data, is implemented in the gstudio R package.

## Key Findings

1. Population Graphs show populations *are* connected, not *which way* gene flow moves — genetic gravity adds direction to that existing, familiar topology, using only the information already in the graph.
2. The core idea: give each population its own local "gravity well," built from the genetic distances to its immediate graph neighbors, then compare how the two populations on either end of an edge each perceive that same edge.
3. The resulting asymmetry index, Δ, is signed: Δ > 0 on edge (i, j) means population *i* is disproportionately shaping the genetic structure of *j* — a source → sink signature; Δ ≈ 0 means balanced, symmetric exchange.
4. No new population-genetic model, no migration–drift equilibrium assumption, no spatial data, and no free parameters — it works from a single snapshot of genotype data.
5. Validated with individual-based simulations of a 25-population stepping-stone metapopulation under three migration regimes: symmetric, and two flavors of directional (flux-conserved and rate-conserved) flow.
6. Most trustworthy as a graph-wide average (Δ̄, and the design-driven axis test) — individual-edge directions are more sensitive to modeling choices (like bandwidth) and are better treated as exploratory leads than firm conclusions.
7. More specific than the existing tool divMigrate: because it reads direction off the graph's conditional structure rather than raw pairwise differentiation, it isn't fooled by genetic drift into reporting "direction" that isn't really there, and it doesn't break down when populations approach fixation.
8. Comes with a matched hierarchy of significance tests — whole-graph, pre-specified axis, and per-edge — plus permutation-based null models for each, all implemented in the [gstudio](https://github.com/dyerlab/gstudio) R package.
