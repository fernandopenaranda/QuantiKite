<img src="logo.png" alt="drawing" width="300"/>

**QuantiKite** is an open-source Julia software which builds a bridge between the two open-source packages of [**Quantica**](https://github.com/pablosanjose/Quantica.jl/tree/master) and [**KITE**](https://github.com/quantum-kite/kite).
It combines the robust customizability of Quantica.jl in constructing quantum tight-binding models with KITE's efficient Chebyshev-inspired algorithms for performing **real-space** bulk spectral and transport calculations in disordered systems. 

In-built multithreading and the on-the-fly nature of KITE routines allow to reach enormous system sizes (up to multi billion atomic orbitals according to KITE's [documentation](https://github.com/quantum-kite/kite)).

KITE uses order-N approximate routines based on the **Kernel Polynomial Method**[^1][^2].

In summary, this package offers a new access to KITE software for the Julia Community.

[^1]: Weiße, A., Wellein, G., Alvermann, A., & Fehske, H. (2006). The kernel polynomial method. [Rev. Mod. Phys., 78, 275–306](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275).
[^2]: M. Joao et al., R. Soc (2020). KITE: high-performance accurate modelling of electronic structure and response functions of large molecules, disordered crystals and heterostructures. [Open Sci. 7, 191809 ](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.78.275).

## Workflow

**KITE** is a compiled (C++) package that takes as input an h5 file containing information about:
  1. Physical properties of the system including its Hamiltonian together with its modifiers (disorder models, (magnetic) vector potential, etc...)
  2. Calculation settings: which observable will be computed and which instructions are needed for the Chebyshev expansion.

This .h5 file is generated through a Python interface that uses [Pybinding](https://docs.pybinding.site/en/stable/) software.

QuantiKite replaces this dependency on Python code by a new API in Julia based on Quantica.jl. It reads as follows: 
The exported function `h5gen` generates a .h5 file taking as inputs: an object `h::Quantica.Hamiltonian` containing all physical information about the periodic part of the system, and julia structs addressing points (1) and (2) above, respectively.

This Hamiltonian object can be regarded as the unit cell that will be repeated in the Bravais lattice directions as many times as determined by the additional settings passed to KITE. Therefore, typically `h` is a small sparse matrix even if the system under study is very large. The non-periodic part of the system is passed directly to KITE using some fields in the settings (see the [Examples](https://github.com/fernandopenaranda/QuantiKite/tree/main/examples) section). 

Once the h5 file is generated it can be passed to KITE's compilers (KITEx and KITE-tools) to perform the desired calculation[^2]. Some plotting functions based on [Makie.jl](https://docs.makie.org/stable/) are also provided for an easy visualization of the KITE outcomes.

## Available functionalities
In this initial version of QuantiKite, we do not have access to the complete catalog of Kite functionalities. At the time of this proof-of-concept release QuantiKite is only suited for  the following calculations in systems with periodic disorder (within the unitcell):

- Density of states (DOS)
- ARPES
- Linear optical conductivities
- Non-linear optical conductivities

Which constitute a comprehensive but not exhaustive list of KITE possibilities.

### Upcoming realeases

A generalization to cover the full functionality of KITE will be soon addressed, including the important feature of non-periodic modifiers (beyond unit-cell disorder) which makes profit of the on-the-fly nature of KITE routines.
   
These include, e.g.: lattice and onsite disorder and Peierls phases as a result of magnetic vector potentials.

## Contributors

[Pablo San-Jose](https://github.com/pablosanjose): Guidance on project design and integration with Quantica.

## Acknowledgements

Funding for this project was obtained through the PID2021-128760NB-I00 project: *Multiscale modeling of twisted bilayer graphene*.

[![Build Status](https://github.com/fernandopenaranda/QuantiKite.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fernandopenaranda/QuantiKite.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fernandopenaranda.github.io/QuantiKite.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fernandopenaranda.github.io/QuantiKite.jl/dev/)
[![Build Status](https://github.com/fernandopenaranda/QuantiKite.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fernandopenaranda/QuantiKite.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/fernandopenaranda/QuantiKite.jl.svg?branch=main)](https://travis-ci.com/fernandopenaranda/QuantiKite.jl)
[![Coverage](https://codecov.io/gh/fernandopenaranda/QuantiKite.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/fernandopenaranda/QuantiKite.jl)
