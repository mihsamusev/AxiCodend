using System;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Globalization;
using CSparse.Storage;
using CSparse;
using CSparse.Double.Factorization;

namespace AxiCodend
{
    class Program
    {
        static void Main(string[] args)
        {
            CultureInfo.DefaultThreadCurrentCulture = new CultureInfo("en-US");

            try
            {
                Console.WindowWidth = 100;
                Console.WindowHeight = 50;
            }
            catch
            {

            }

            //=========================
            // INPUT
            //=========================

            //string path = GetExePath();
            //string input = path + "\\input.txt";

            var path = new PathsIO();

            if (!File.Exists(path.input))
            {
                Console.WriteLine("input.txt does not exis in the assembly path.");
                return;
            }

            var hexMeshPanel = new HexMeshPanelMaterial(path);
            hexMeshPanel.PrintInfo();

            AxiCodend codend = new AxiCodend();
       
            if (hexMeshPanel.MeshOrientation == 0)
            {
                codend = new AxiModelT0(path);   // To start with initialize T0 model
            }
            else if(hexMeshPanel.MeshOrientation == 90)
            {
                //codend = new AxiModelT90(path);   // To start with initialize T0 model
            }
            else
            {
                throw new IOException("Incorrect mesh orientation. Value should be 0 or 90.");
            }
            
            var towing = new Towing(path);

            var catches = new Catch(path);

            var TowingSimulation = new Simulation(codend, catches, towing)
            {
                SolverSettings = new SolverSettings(path),
                Paths = path
            };

            TowingSimulation.Simulate();

            //if (!Console.IsOutputRedirected)
            //{
            //    Console.WriteLine("Press any key to close.");
            //    Console.ReadKey();
            //}
        }
        
        // Static methods 
        //=======================================================================

        static string GetExePath()
        {
            string fullPath = Path.GetDirectoryName(Assembly.GetExecutingAssembly().GetName().CodeBase);
            return fullPath.Replace("file:\\", "");
        }
    
        static void SaveArrayToPath(double[] X, string path)
        {
            using (StreamWriter ResultFile = new StreamWriter(path))
            {
                for (int i = 0; i < X.Length; i++)
                {
                    ResultFile.WriteLine(X[i]);
                }
            }
        }

        static void SaveResultsToPath(string[] Results, string path)
        {

        }

        static void IterationResults(int iter,int restarts, Stopwatch stopwatch)
            {
                Console.WriteLine("\nConvergence acheived in {0} " + "" +
                                    "iterations and {1} restarts.\n" +
                                    "Total time elapsed: {2:mm\\:ss} [min:sec]\n",
                                    iter, restarts, stopwatch.Elapsed);
            }

        static public double Parab3p(double lc, double lp, double fc, double fp, double f0)
        {
            /* Use line search to find out lambda that decreases the residual R

            X = X + lambda*h

            If full step (lambda = 1) and half step (lambda = 0.5) are rejected 
            due to insufficiend decrease in comparison with the residual from previous NR iteration
            next line search iterations calculate lambda using 3-point parabolic fitting with safeguard
            In addition, if pre-last iteration is neither succesful, the final step is chosen according 
            to step limit method (either full step or a fraction of a twine length depending on size of vector h)
            */

            /*      input:
            lc = current steplength
            lp = previous steplength

            fc = value of || F(x_c + lc d) ||^ 2
            fp = value of || F(x_c + lp d) ||^ 2
            f0 = value of || F(x_c)        ||^ 2
             */
            const double sigma0 = 0.1;
            const double sigma1 = 0.5;
            /*
            Compute coefficients of interpolation polynomial
            p(l) = ff0 + (c1 l + c2 l ^ 2) / d1
            p(l)' = (c1 + 2 c2 l) / d1
            p(l)'' = 2 c2 / d1

            where:
            c1 = lc*lc*(fp-f0) - lp*lp*(fc-f0);
            c2 = lp*(fc-f0) - lc*(fp-f0);
            d1 = (lc - lp) * lc * lp < 0

            so if c2 > 0 then p(l)'' < 0 we have negative curvature
            and default to lp = sigma1 * l
            */

            double c2 = lp * (fc - f0) - lc * (fp - f0);

            if (c2 >= 0)
            {
                return sigma1 * lc;
            }

            double c1 = lc * lc * (fp - f0) - lp * lp * (fc - f0);
            double lmin = - c1 * 0.5 / c2;

            /* safeguard */

            if (lmin < sigma0 * lc)
            {
                lmin = sigma0 * lc;
            }

            if (lmin > sigma1 * lc)
            {
                lmin = sigma1 * lc;
            }

            return lmin;
        }

        static void SparseSolveSystem(CompressedColumnStorage<double> Ksparse, double[] F, ref double[] h)
        {
            var lu = SparseLU.Create(Ksparse, ColumnOrdering.MinimumDegreeAtPlusA, 1);
            lu.Solve(F, h);
        }

        static void DisplaySparseMatrix(CompressedColumnStorage<double> S)
        {
            Console.WriteLine();
            foreach (var Triple in S.EnumerateIndexed())
            {
                Console.WriteLine("({0}, {1})\t{2}", Triple.Item1, Triple.Item2, Triple.Item3);
            }
            Console.WriteLine();
        }

        static void DisplayFullMatrix(CompressedColumnStorage<double> A, int iStart, int jStart, int iEnd, int jEnd)
        {
            for (int i = iStart; i <= iEnd; i++)
            {
                for (int j = jStart; j <= jEnd; j++)
                {
                    Console.Write("{0,20:F4}", A.At(i, j));
                }
                Console.Write("\n");
            }
        }

        static void DisplayVector(double[] A)
        {
                for (int j = 0; j < A.Length; j++)
                {
                    Console.WriteLine("{0,10:F4}", A[j]);
                }
        }
    }
}