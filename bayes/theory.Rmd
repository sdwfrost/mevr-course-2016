# Bayesian inference of time trees

## Introduction

- Bayesian inference tries to obtain a probability distribution rather than a single best fit c.f. maximum likelihood
- This is typically obtained using numerical approaches such as Markov Chain Monte Carlo
- For example, Metropolis sampling involves the following
  - Make a change to the parameters
  - If the new parameters result in a better fit, change to them
  - If the new parameters result in a worse fit, move to them with a probability based on the drop in fit
  - Repeat many times

## Strict and molecular clocks

- In phylogenetic inference, Bayesian approaches have become popular to fit molecular clocks to sequence data
- Two commonly used platforms
  - MrBayes
  - BEAST
- Both use the same underlying principles, but the models are parameterised a little differently
  - BEAST has a wide range of 'demographic models' that can be used as prior distributions for the branch lengths
- The input is also different
  - BEAST has an XML input format
  - MrBayes uses the Nexus format
