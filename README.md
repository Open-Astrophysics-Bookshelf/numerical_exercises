# Computational Hydrodynamics for Astrophysics

*part of the Open Astrophysics Bookshelf*

Notes on numerical methods for computational astrophysical hydrodynamics.

These notes describe the way I think about the numerical methods commonly
used with grid-based codes in astrophysical hydrodynamics.  The notes
are written in LaTeX, and should build by typing 'make' in the main
directory.

Working implementations for all of the solvers are contained either in
the main pyro code or in the set of hydro examples, both referenced
below.

## Chapters

The following chapters are mostly written:

- Simulation Overview
- Finite-Volume Grids
- Advection
- Burgers' Equation
- Euler Equations: Theory
- Euler Equations: Numerical Methods
- Elliptic Equations and Multigrid
- Diffusion
- Model Multiphysics Problems
- Reactive Flow
- Planning a Simulation
- Incompressible Flow and Projection Methods
- Low Mach Number Methods

The following are things I'd like to add in the next 1-2 years:

- Fluid Instabilities
- Rotation and Self-gravity
- Radiation Hydrodynamics
- MHD
- AMR
- Mapped Grids

The following are things hopefully will eventually get written:

- Relativisitc Flows
- Higher-Order Methods
- Implicit Hydrodynamics


## Getting PDFs and Source

A PDF version of these notes is available at:

http://bender.astro.sunysb.edu/hydro_by_example/CompHydroTutorial.pdf


There are two sets of companion codes that go along with these notes:

 - *hydro_examples*: https://github.com/zingale/hydro_examples

   simple, standalone, 1-d solvers that illustrate the basic ideas.

 - *pyro*: https://github.com/zingale/pyro2

   pyro is a 2-d full simulation code designed with simplicity in
   mind that implements the core solvers described in these notes
   along with various test problems.




