using System.Diagnostics;
using CSparse.Double.Factorization;
using CSparse;

namespace AxiCodend
{
    class Simulation
    {
        //============================
        // VARIABLES
        //============================

        public AxiCodend Codend;
        public double towing_speed;
        public int[] catches;


        public SolverSettings SolverSettings { get; set; }
        public ICodendSaver ResultSaver {get; set;}

        private double[] precalcX;
        private double[] previousX;

        //============================
        // CONSTRUCTOR
        //============================

        public Simulation(AxiCodend Codend, int[] catches, double towing_speed, ICodendSaver result_saver)
        {
            this.Codend = Codend;
            this.catches = catches;
            this.towing_speed = towing_speed;

            this.ResultSaver = result_saver;
            SolverSettings = new SolverSettings(); // go with default settingst
            precalcX = new double[Codend.dof];
            previousX = new double[Codend.dof];
        }

        //============================
        // METHODS
        //============================

        private double DisplacmentNorm(double[] h)
        {
            double norm = 0;
            for (int i = 0; i < Codend.dof; i++)
            {
                norm += Math.Pow(h[i], 2);                      // residue as sum of squares
            }
            return Math.Sqrt(norm);
        }

        private bool CheckPrecalc()
        {
            /*Check if precalculation of the initial shape is requred*/
            bool usePrecalc = false;

            if (Codend.GetMeshesBlocked() < SolverSettings.MinCatchBlock) // due to low catch
            {
                Codend.ApplyCatch(SolverSettings.MinCatchBlock);
                Console.WriteLine("\nCodend is too empty, using {0} meshes blocked by catch instead to precalculate initial shape.\n",SolverSettings.MinCatchBlock);
                usePrecalc = true;
            }

            if (Codend.towSpeed < SolverSettings.MinTowingSpeed)                      // due to very slow towing
            {
                Codend.ApplyTowing(SolverSettings.MinTowingSpeed);
                Console.WriteLine("\nPrecalculating initial shape using speed of {0,4:F3} m/s...\n",SolverSettings.MinTowingSpeed);
                usePrecalc = true;             
            }

            return usePrecalc;
        }

        private void ReportRestartCause(bool hasDiverged)
        {
            if (hasDiverged)
            {
                Console.WriteLine("\nIteration diverges or the algebraic system is too ill-conditioned." +
                                  "\nRestarting with higher added diagonal stiffness.\n");
            }
            else
            {
                Console.WriteLine("\nConverged to incorrect tangled solution." +
                                  "\nRestarting with higher added diagonal stiffness.\n");
            }
        }

        private void SetStageInitialShape(int stage)
        {
            if (stage == 2)
            {
                precalcX.CopyTo(Codend.X, 0);
                return;
            }

            if (SolverSettings.UsePreviousAsPrecalc && DisplacmentNorm(previousX) != 0)
            {
                previousX.CopyTo(Codend.X, 0);
            }
            else
            {
                Codend.SetInitialShapeSmooth(); 
            }
        }

        private bool IsIncorrectSolution()
        {
            if (Codend.X.Min() < 0)
            {
                return true;
            }

            double Ranalytical = 2 * Codend.GetMeshesAround() * 0.5 * Codend.Material.MeshSide / (Math.PI * Math.Sqrt(6));
            double tol = 2;

            if (Codend.MaxRadius() > tol * Ranalytical)
            {
                return true;
            }

            double stretchedLength = Codend.GetMeshesAlong() * (Codend.Material.MeshSide + 2 * Codend.Material.KnotSize);
            if (Codend.Length() > stretchedLength)
            {
                return true;
            }

            return false;
        }

        private void Solve()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            double[] h = new double[Codend.dof];
            double PrevTowSpeed = Codend.towSpeed;
            int PrevCatch = Codend.GetMeshesBlocked();

            SparseLU LU;
            var LUorder = ColumnOrdering.MinimumDegreeAtPlusA;
            double LUtolerance = 1.0;

            bool hasConverged = false;                          // convergence flag
            bool correctSolution = false;                       // correctness of solution flag
            bool hasDiverged = false;                           // restart flag
            bool usePrecalc = CheckPrecalc();
            int restartsCount = 0;
            int iter = 0;

            double hNorm = 1;
            double addStiff = SolverSettings.DiagStiffness;
            double addStiff0 = SolverSettings.DiagStiffness;
            double tolStiff = SolverSettings.StiffnessTol;

            //=============================         
            // Start calculation
            //=============================

            for (int stage = 1; stage <= 2; stage++)
            {
                //===========================
                // Restarts loop
                //===========================

                while (!correctSolution)
                {
                    // Set initial shape depending on whether it is precalculated from previous stage or not
                    Codend.ClearState();

                    SetStageInitialShape(stage);

                    Codend.UpdateTotalForces(IncludeBC: true);            
                    Codend.UpdateResidual();

                    if (addStiff0 == 0 && restartsCount > 0)  // if the scheme doesnt work without extra stiffness  
                    {
                        addStiff0 = 1;
                    }
                    addStiff = addStiff0;
                    tolStiff = SolverSettings.StiffnessTol;
                    //===========================
                    // Newton-Raphson loop
                    //===========================

                    while (!hasConverged)
                    {
                        // check if diverged
                        if (Codend.R > SolverSettings.ResidualMax)
                        {
                            hasDiverged = true;
                            break;
                        }

                        // check whether to add less diagonal stiffness
                        if (Codend.R < tolStiff)               
                        {
                            addStiff *= SolverSettings.ReduceStiffnessBy;
                            tolStiff *= SolverSettings.ReduceStiffnessBy;
                        }

                        // check if converged and no stiffness is added
                        if (Codend.R < SolverSettings.ResidualTol && hNorm < SolverSettings.DisplacementTol)
                        {
                            break;
                        }
                       
                        Codend.UpdateTotalJacobian(kDiag: -addStiff, IncludeBC: true);

                        LU = SparseLU.Create(Codend.J, LUorder, LUtolerance);
                        LU.Solve(Codend.F, h);

                        hNorm = DisplacmentNorm(h);
                        Codend.UpdatePosition(h, 1);

                        Codend.UpdateTotalForces(IncludeBC: true);
                        Codend.UpdateResidual();

                        iter++;                                                   // display iteration status
                        Console.WriteLine(
                            PrintIter(iter, 0, addStiff, Codend.R, hNorm));
                    }

                    //===================
                    // Restart routine               
                    //===================

                    if (hasDiverged || IsIncorrectSolution())
                    {
                        restartsCount += 1;
                        addStiff0 *= SolverSettings.IncreaseStiffnessBy;
                        ReportRestartCause(hasDiverged);
                        hasDiverged = false;
                    }
                    else
                    {
                        correctSolution = true; // correct solution found
                        Console.WriteLine("\nConvergence acheived in {0} iterations and {1} restarts in {2} [ms]\n",
                        iter, restartsCount, stopwatch.ElapsedMilliseconds);
                    }
                }


                if (usePrecalc)
                {
                    Codend.X.CopyTo(precalcX, 0);
                    Codend.ApplyCatch(PrevCatch);
                    Codend.ApplyTowing(PrevTowSpeed);

                    Console.WriteLine("\nInitial shape is precalculated. Moiving to and actual case.\n");
                    correctSolution = false;
                    usePrecalc = false;
                    hasConverged = false;
                }
                else
                {
                    break;
                }
            }
        }

        public void Simulate()
        {
            Stopwatch stopwatch = new Stopwatch(); // can wrap entire Solver class
            stopwatch.Start();

            double[][] allX = new double[catches.Length][];
            ResultSaver.AddHeader(Codend.Material, Codend.Geometry);

            for (int i = 0; i < catches.Length; i++)
            {
                Console.WriteLine(
                    "Simulation {0} with {1} blocked meshes\n", i, catches[i]);

                Codend.ApplyTowing(towing_speed);
                Codend.ApplyCatch(catches[i]);
                Solve();

                ResultSaver.AddRun(
                    catches[i],
                    towing_speed,
                    Codend.GetMetrics(),
                    Codend.X);

                if (SolverSettings.UsePreviousAsPrecalc)
                {
                    Codend.X.CopyTo(previousX, 0);
                }

                allX[i] = Codend.X;                
            }


            Console.WriteLine(
                "\nSimulation of {0} catches is finished in {1} [ms]\n",
                catches.Length,
                stopwatch.ElapsedMilliseconds);

            ResultSaver.Finish();

        }

        private string PrintIter(int iter, int funEval, double addStiff, double R, double hNorm)
        {
            return String.Format("Iter: {0,-7:D}" +
                    "LS steps: {1,-6:D}" +
                    "addStiff: {2,-10:e2}{3,-8:S}" +
                    "R: {4,-10:e2}{5,-6:S}" +
                    "|h|: {6,-10:e2}{7,-6:S}", iter, funEval, addStiff, "N/m", R,"N", hNorm,"m");
        }
    }
}
