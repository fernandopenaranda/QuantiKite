# QuantiKite: A wrapper between [Quantica](https://github.com/pablosanjose/Quantica) and [Kite](https://github.com/quantum-kite/kite)

## Description

QuantiKite aims to build a bridge between the two wonderful open-source packages of Quantica.jl[^1] and KITE[^2].
QuantiKite.jl combines the robust customizability of Quantica.jl in constructing Quantum tight-binding models with efficient Chebyshev-inspired algorithms in KITE for performing real-space bulk spectral calculations in disordered systems. In-built multithreading and the on-the-fly nature of KITE routines allows to reach enormous system sizes up to multi billion atomic orbitals[^2].

Kite routines are based on order N expansions known as the Kernel Polynomial Method (KPM)[^3] and its documentation can be found in [. M. Joao et al., R. Soc. Open Sci. 7, 191809 (2020) [https://royalsocietypublishing.org/doi/full/10.1098/rsos.191809].

This package is purely written in Julia as Quantica, so there is no need to interface with different programming languages.

## Workflow

Kite is a compiled (c++) package that takes as input for its calculations a h5 file containing information about
  1. Physical properties of the system including its Hamiltonian together with its modifiers (disorder models, (magnetic) vector potential, etc
  2. Settings for the precise type of calculation (dos, optical conductivities) and instructions for the Chebyshev expansion

QuantiKite aims to write under a friendly API this h5 file required by KITE using information coming from Quantica. More precisely, it provides a function `h5gen` which takes an object `h::Quantica.Hamiltonian` containing all physical information about the periodic part of the system and julia structs addressing the points (1) and (2) above, respectively.

This Hamiltonian object contains the information about the periodic part of the system and can be regarded as the unit cell that will be repeated in the Bravais lattice directions as many times as determined by the additional settings passed to Kite (see [] or []). Therefore, typically h is a small sparse matrix even if the system under study is very large. The non-periodic part of the system is passed directly to Kite using some fields in the settings (see Examples in or in for more information). 

Once the h5 file is generated it can be passed to Kitex and Kite-tools compiled functions to get the result. Some plotting functions based on Makie.jl[^4] are also provided under this package for an easy visualization of the Kite outcomes.

## Available functionalities
In this first version of QuantiKite we cannot access to the whole catalogue of Kite functionalities. At the moment of this proof-of-concept release we can already compute in periodic systems:

- Density of states (DOS)
- Linear optical conductivities
- Non-linear optical conductivities

  In the current release we can already obtain these observables in non-periodic disordered systems or systems subjected to magnetic fields using Kite. However, this is done by enlarging the unticell to cover the whole system and introducing the disorder through Quantica in h and at the cost of not accessing the on-the-fly routines in KITE.
  This can be regarded as if the supercell equals the unitcell and no modifiers are required. This capability will be soon corrected.
   
#### Soon

- Non-periodic modifiers including (random) lattice, onsite, and hopping disorder or Peierls phases as a result of magnetic vector potentials.
- See the complete observable list in KITE docs

## Acknowledgments

##

[![Build Status](https://github.com/fernandopenaranda/QuantiKite.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fernandopenaranda/QuantiKite.jl/actions/workflows/CI.yml?query=branch%3Amain)

https://github.com/fernandopenaranda/ChebyshevExpansions.git


# QuantiKite

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fernandopenaranda.github.io/QuantiKite.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fernandopenaranda.github.io/QuantiKite.jl/dev/)
[![Build Status](https://github.com/fernandopenaranda/QuantiKite.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fernandopenaranda/QuantiKite.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/fernandopenaranda/QuantiKite.jl.svg?branch=main)](https://travis-ci.com/fernandopenaranda/QuantiKite.jl)
[![Coverage](https://codecov.io/gh/fernandopenaranda/QuantiKite.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/fernandopenaranda/QuantiKite.jl)
