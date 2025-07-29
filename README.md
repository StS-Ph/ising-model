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

## The one-dimensional Ising Model
The IMTF describes the nearest neighbor spin coupling on a linear chain, biased by an external transverse magnetic field. It experiences a quantum phase transition based on the relation of the spin coupling *J* and the magnetic field parameter *h*, which can be studied using this simulation. A detailed discussion, including the implemented Hamiltonian, can be found in *P. Pfeuty. Ann. Phys. 57 1.* and [here](https://en.wikipedia.org/wiki/Transverse-field_Ising_model).

## Lanczos Algorithm
The Lanczos Algorithm can be used to reduce complexity of diagonalizing a full Hamiltonian matrix (dimensions scale exponentially in the system size). 
It projects the Hamiltonian onto a Krylov subspace, reducing the dimension down to the algorithm iterations. The resulting matrix is additionally of tri-diagonal form
and can be very easily diagonalized. A detailed discussion can be found [here](https://www.cond-mat.de/events/correl15/manuscripts/koch.pdf).
