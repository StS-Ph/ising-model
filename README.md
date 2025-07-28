# Numerical simulation of the one-dimensional Ising model in a transversal field using the Lanczos algorithm

This repository contains the C# code for the implementation of the Lanczos algorithm and the one-dimensional Ising model in a transversal field (IMTF). 
Additionally, relevant classes and methods for the diagonalization of the IMTF in different parameter regimes using the Lanczos algorithm are implemented.

## Getting started
The class `SimulationRunner` provides the static method `RunParameterVariation`, which can be used to diagonalize the Ising Hamiltonian in different parameter regimes.
Different coupling parameters can here be submitted as a vector. The simulation results are thereby directly saved to `.txt` files.

The class `IsingHamiltonian` implements the IMTF and the Hamiltonian matrix can be generated via the `BuildHamiltonian` method. 
The Lanczos algorithm is implemented with the static class `LanczosSolver`, 
simply providing the Lanczos algorithm followed by remaining diagonlization via the static method `ComputeEigenvalues`.
The saving of simulation results to file is handled by the `FileExporter` class.
