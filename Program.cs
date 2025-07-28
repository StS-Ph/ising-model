/*
This program implements the Lanczos Algorithm, a numerical method to convert a Hamiltonian into a tri-diagonal form.
This tri-diagonal Hamiltonian can then be easier diagonalized compared to the original Hamiltonian.
Additionally, the Hamiltonian of the one-dimensional Ising Modell in a transversal magnetic field is implemented.
The numerical simulations apply the Lanczos algorithm to this Hamiltonian for different parameter regimes.
*/
using System;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;

namespace Ising_App
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Start the calculations by typing enter.");
            Console.ReadLine();

            int lanIt = 30;
            int n = 10;
            double[] qhJValues = { 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 };

            SimulationRunner runner = new SimulationRunner();
            runner.RunParameterVariation(qhJValues, n, lanIt);

            Console.WriteLine("Finished.");
            Console.ReadLine();
        }
    }

    class SimulationRunner
    {
        public void RunParameterVariation(double[] couplings, int systemSize, int lanczosIterations)
        {
            double[] e0Vec = new double[couplings.Length];
            double[] e1Vec = new double[couplings.Length];
            double[] gapVec = new double[couplings.Length];

            // Running the Lanczos algorithm for a fixed system size and a fixewd number of iterations but with varying coupling strength
            for (int i = 0; i < couplings.Length; i++)
            {
                double qhJ = couplings[i];

                // initialize Ising Hamiltonian with current coupling strength
                IsingHamiltonian hamiltonian = new IsingHamiltonian(systemSize, qhJ);
                var matrix = hamiltonian.BuildHamiltonian();
                // initial state to find ground state
                var groundState = SpinState.FromIndices(systemSize, new int[] { });
                // initial state to find first excitation
                var excitedState = SpinState.FromIndices(systemSize, new int[] { 0 });

                // Lanczos Algorithm for ground state
                var eigVals = LanczosSolver.ComputeEigenvalues(matrix, groundState, lanczosIterations);
                // ground state energy
                e0Vec[i] = eigVals[0];

                // Lanczos Algorithm for first excited state
                eigVals = LanczosSolver.ComputeEigenvalues(matrix, excitedState, lanczosIterations);
                // energy
                e1Vec[i] = eigVals[0];

                // Calculate energy gap
                gapVec[i] = e1Vec[i] - e0Vec[i];

                Console.WriteLine($"QhJ={qhJ}, E0={e0Vec[i]}, E1={e1Vec[i]}, gap={gapVec[i]}");
            }

            FileExporter.WriteToTxt("E0_vs_QhJ", couplings, e0Vec, "QhJ", "E0");
            FileExporter.WriteToTxt("E1_vs_QhJ", couplings, e1Vec, "QhJ", "E1");
            FileExporter.WriteToTxt("Gap_vs_QhJ", couplings, gapVec, "QhJ", "E1-E0");
        }
    }
    static class Pauli
    {
        // Provide the relevant Pauli matrices
        private static readonly Matrix<double> sigmaX = Matrix<double>.Build.DenseOfArray(new double[,]
        {
            { 0, 1 },
            { 1, 0 }
        });

        private static readonly Matrix<double> sigmaZ = Matrix<double>.Build.DenseOfArray(new double[,]
        {
            { 1, 0 },
            { 0, -1 }
        });

        public static Matrix<double> X => sigmaX;
        public static Matrix<double> Z => sigmaZ;
    }
    class IsingHamiltonian
    {
        // Class to represent the Ising Hamiltonian in the transversal magnetic field
        public int SystemSize { get; } // number of spins
        public double Coupling { get; } // coupling strength between spins and magnetic field

        public IsingHamiltonian(int n, double c)
        {
            SystemSize = n;
            Coupling = c;
        }
        public Matrix<double> BuildHamiltonian()
        {
            // Create Hamiltonian matrix of the Ising Hamiltonian in the transversal magnetic field with set number of spins and coupling strength

            // Get relevant Pauli matrices
            Matrix<double> PauliZ = Pauli.Z;
            Matrix<double> PauliX = Pauli.X;
            // create operator for two-body XX interaction of neighbouring spins
            Matrix<double> OpXX = PauliX.KroneckerProduct(PauliX);

            var M = Matrix<double>.Build;
            int dim = Convert.ToInt32(Math.Pow(2, SystemSize));
            Matrix<double> sum1 = M.SparseDiagonal(dim, dim, 0.0);
            Matrix<double> sum2 = M.SparseOfMatrix(sum1);

            // Create magnetic field term
            for (int i = 0; i <= SystemSize - 1; i++)
            {
                int lExp = Convert.ToInt32(Math.Pow(2, i));
                int rExp = Convert.ToInt32(Math.Pow(2, SystemSize - 1 - i));
                Matrix<double> x = M.SparseIdentity(lExp).KroneckerProduct(PauliZ);
                Matrix<double> y = x.KroneckerProduct(M.SparseIdentity(rExp));
                sum1 = sum1.Add(y);
            }

            // Create two-body interaction term
            for (int i = 0; i <= SystemSize - 2; i++)
            {
                int lExp = Convert.ToInt32(Math.Pow(2, i));
                int rExp = Convert.ToInt32(Math.Pow(2, SystemSize - 2 - i));
                Matrix<double> x = M.SparseIdentity(lExp).KroneckerProduct(OpXX);
                Matrix<double> y = x.KroneckerProduct(M.SparseIdentity(rExp));
                sum2 = sum2.Add(y);
            }

            /* cyclic chain (periodic boundary conditions */
            int cycExp = Convert.ToInt32(Math.Pow(2, SystemSize - 2));
            Matrix<double> cyc1 = PauliX.KroneckerProduct(M.SparseIdentity(cycExp));
            Matrix<double> cyc2 = cyc1.KroneckerProduct(PauliX);
            sum2 = sum2.Add(cyc2);

            // add coupling factor
            sum1 = sum1.Multiply(Coupling);
            Matrix<double> H = sum1.Add(sum2);
            H = H.Multiply(-1.0);

            return H;

        }
    }
    static class LanczosSolver
    {
        // Class to represent the Lanczos Algorithm
        public static Vector<double> ComputeEigenvalues(Matrix<double> H, Vector<double> initialState, int iterations)
        {
            // Applying the Lanczos algorithm to a given Hamiltonian H with an certain initialState
            int n = iterations; // number of algorithm iterations
            double[] dArr = new double[n];
            double[] fArr = new double[n];
            Vector<double> newState = initialState.Clone();
            Vector<double> currState = initialState.Clone();
            Vector<double> oldState = initialState.Clone();

            Vector<double> gamma = H.Multiply(currState);
            dArr[0] = currState.DotProduct(gamma);
            gamma = gamma.Subtract(currState.Multiply(dArr[0]));

            fArr[0] = gamma.DotProduct(gamma);
            fArr[0] = Math.Sqrt(fArr[0]);
            newState = gamma.Divide(fArr[0]);

            currState = newState;
            //double test;
            //double test2;
            for (int i = 1; i <= n - 1; i++)
            {
                gamma = H.Multiply(currState);
                //test = gamma.DotProduct(gamma);
                //test2 = currState.DotProduct(currState);
                dArr[i] = currState.DotProduct(gamma);
                gamma = gamma.Subtract(currState.Multiply(dArr[i]));
                //test = gamma.DotProduct(gamma);
                gamma = gamma.Subtract(oldState.Multiply(fArr[i - 1]));
                //test = gamma.DotProduct(gamma);

                fArr[i] = gamma.DotProduct(gamma);
                fArr[i] = Math.Sqrt(fArr[i]);

                newState = gamma.Divide(fArr[i]);

                oldState = currState;
                currState = newState;
            }

            Matrix<double> HTriDiag = Matrix<double>.Build.SparseOfDiagonalArray(dArr);
            HTriDiag[1, 0] = fArr[0];
            for (int i = 1; i <= n - 2; i++)
            {
                HTriDiag[i - 1, i] = fArr[i - 1];
                HTriDiag[i + 1, i] = fArr[i];
            }
            HTriDiag[n - 2, n - 1] = fArr[n - 2]; // Final tridiagonal matrix from Lanczos 

            Evd<double> evd = HTriDiag.Evd(); // use ed to diagonalize it
            // return sorted vector of eigenvalues of tridiagonal Hamiltonian
            var evals = evd.EigenValues.Real().ToArray();
            Array.Sort(evals);
            return Vector<double>.Build.DenseOfArray(evals);
        }
    }
    class SpinState
    {
        // Class to represent spin states as vectors
        public static Vector<double> FromIndices(int N, int[] indices)
        {
            //This Method calculates the spin state vector, using the kronecker product.
            //Make sure counting starts at zero, not one for indices!
            if (indices.Length == 0)
            {
                //Ground state
                int dim = Convert.ToInt32(Math.Pow(2, N));
                Vector<double> stateVec = Vector<double>.Build.Sparse(dim, 0.0);
                stateVec[0] = 1.0;
                return stateVec;
            }
            else if (indices.Length == N)
            {
                // fully excited state
                int dim = Convert.ToInt32(Math.Pow(2, N));
                Vector<double> stateVec = Vector<double>.Build.Sparse(dim, 0.0);
                stateVec[dim - 1] = 1.0;
                return stateVec;
            }
            var M = Matrix<double>.Build;
            Matrix<double> spinDown = M.Sparse(2, 1, 0.0);
            Matrix<double> spinUp = M.Sparse(2, 1);
            spinDown[0, 0] = 1.0;
            spinUp[1, 0] = 1.0;

            Matrix<double> spinStateM = M.SparseOfMatrix(spinDown);
            int cnt = 0;
            if (indices[cnt] == 0)
            {
                spinStateM = M.SparseOfMatrix(spinUp);
                cnt++;
            }

            for (int i = 1; i <= N - 1; i++)
            {
                if (cnt >= indices.Length)
                {
                    spinStateM = spinStateM.KroneckerProduct(spinDown);
                }
                else if (i == indices[cnt])
                {
                    spinStateM = spinStateM.KroneckerProduct(spinUp);
                    cnt++;
                }
                else
                {
                    spinStateM = spinStateM.KroneckerProduct(spinDown);
                }
            }

            double[] spinStateArr = spinStateM.ToColumnMajorArray();
            Vector<double> spinStateVec = Vector<double>.Build.SparseOfArray(spinStateArr);

            return spinStateVec;
        }
    }
    static class FileExporter
    {
        // Class for exporting simulation data to file
        public static void WriteToTxt(
            string title,
            double[] xValues,
            double[] yValues,
            string xLabel,
            string yLabel,
            string xUnits = "",
            string yUnits = "",
            string path = "")
        {
            // Method to print all given data in a txt file

            // setup correct filename
            string filename = "";
            if (path == "")
            {
                filename = string.Concat(title, ".txt");
            }
            else
            {
                filename = string.Concat(path, title, ".txt");
            }

            // Write data to file
            using (var writer = new StreamWriter(filename, false))
            {
                writer.WriteLine($"{xLabel} in [{xUnits}];{yLabel} in [{yUnits}]");

                for (int i = 0; i < xValues.Length; i++)
                {
                    writer.WriteLine($"{xValues[i]};{yValues[i]}");
                }
            }
        }
    }
    
} 