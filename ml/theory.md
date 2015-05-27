# Theory

## Introduction

- Methods that fit a model of molecular evolution to the sequence data are more computationally intensive, but typically show better performance than distance based methods

## How many trees?

---

- How many trees are possible with 12 sequences
- [ ] 105
- [ ] 10,295
- [ ] 2,027,025
- [x] 654,729,025

---

## Maximum likelihood inference

- The large number of trees makes it difficult to find the tree with the highest likelihood
- Phylogeny programs have to use heuristic approaches to find the 'best' tree
  - Starting with an initial tree, make modifications and test whether they give a better tree or not
    - Nearest neighbour interchange
    - Subtree pruning and regrafting

## Models of DNA evolution

- Just like with distance based methods, we assume a model of sequence evolution
- However, we can now test which is the 'best' model
- We often don't need a good tree for this

## Rate heterogeneity

- In addition to assuming a model of how, for example, one nucleotide changes to another, we can also assume a model of how substitution rates overall vary across the sequence
  - Constrained regions e.g. those that are functionally important
  - Variable regions e.g. those under immune selection
- Two sorts of models of rate heterogeneity
  - Gamma distribution (possibly with an additional invariant category)
  - Categorical

## Balancing model fit and complexity

- To choose a model, we have to balance model fit (likelihood) with complexity (number of parameters)
- For non-nested models, two criteria are commonly used (**lower is better**)
  - [Akaike's Information Criterion](http://en.wikipedia.org/wiki/Akaike_information_criterion) (AIC)
    - $AIC = 2 k - 2 \ln (L)$
  - [Bayesian Information Criterion](http://en.wikipedia.org/wiki/Bayesian_information_criterion) (BIC)
    - $BIC = k \ln (n) - 2 \ln (L)$
- BIC favours simpler models than AIC
