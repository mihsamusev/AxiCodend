using System;
using System.IO;

namespace AxiCodend
{
    class SolverSettings
    {
        //============================
        // fields
        //============================

        public int IterMax { get; set; }
        public double ResidualTol { get; set; }                   
        public double DisplacementTol { get; set; }                     
        public double ResidualMax { get; set; }              // Maximum residual after which the scheme is considered divergent
        public double DiagStiffness { get; set; }
        public double StiffnessTol { get; set; }
        public double ReduceStiffnessBy { get; set; }
        public double IncreaseStiffnessBy { get; set; }
        public bool ShowLineSearchSteps { get; set; }
        public int LineSearchMax { get; set; }              // maximum amount of inexact line search iterations
        public double AlphaMax { get; set; }                // alpha for Armijo step limit globalization method (20% of twine length)  
        public int MinCatchBlock { get; set; }
        public double MinTowingSpeed { get; set; }
        public bool UsePreviousAsPrecalc { get; set; }

        //============================
        // constructor
        //============================

        public SolverSettings()
        {
            LoadDefault();
        }

        public SolverSettings(PathsIO path)
        {
            LoadDefault();
            LoadInput(path.input);
        }

        //============================
        // methods
        //============================


        private void LoadDefault()
        {
            IterMax = 5000;
            ResidualTol = 1e-3;                          // norm force < 0.001 N
            DisplacementTol = 1e-4;                      // norm dispacement < 0.1 mm
            ResidualMax = Math.Pow(10, 20);              // Maximum residual after which the scheme is considered divergent
            StiffnessTol = 1;            
            DiagStiffness = 1;
            
            ReduceStiffnessBy = 0.1;
            IncreaseStiffnessBy = 2;

            ShowLineSearchSteps = false;
            LineSearchMax = 6;                           // maximum amount of inexact line search iterations
            AlphaMax = 0.2 * 0.16;                       // alpha for Armijo step limit globalization method (20% of twine length)  

            MinCatchBlock = 5;
            MinTowingSpeed = 0.1;
            UsePreviousAsPrecalc = false;
        }

        private void LoadInput(string inpPath)
        {
            string[] names = { "IterMax", "ResidualTol", "DisplacementTol", "ResidualMax", "StiffnessTol",
                               "DiagStiffness", "ReduceStiffnessBy", "IncreaseStiffnessBy", "ShowLineSearchSteps",
                               "LineSeachMax","AlphaMax","MinCatchBlock","MinTowingSpeed","UsePreviousAsPrecalc"};
            string[] lines = File.ReadAllLines(inpPath);
            string[] parts;
            int currentLine = 0;

            foreach (var line in lines)
            {
                parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);

                if (line.Contains(names[0]))
                    IterMax = Convert.ToInt32(parts[1]);

                if (line.Contains(names[1]))
                    ResidualTol = Convert.ToDouble(parts[1]);

                if (line.Contains(names[2]))
                    DisplacementTol = Convert.ToDouble(parts[1]);

                if (line.Contains(names[3]))
                    ResidualMax = Convert.ToDouble(parts[1]);

                if (line.Contains(names[4]))
                    StiffnessTol = Convert.ToDouble(parts[1]);

                if (line.Contains(names[5]))
                    DiagStiffness = Convert.ToDouble(parts[1]);

                if (line.Contains(names[6]))
                    ReduceStiffnessBy = Convert.ToDouble(parts[1]);

                if (line.Contains(names[7]))
                    IncreaseStiffnessBy = Convert.ToDouble(parts[1]);

                if (line.Contains(names[8]))
                    ShowLineSearchSteps = Convert.ToBoolean(Convert.ToInt16(parts[1]));

                if (line.Contains(names[9]))
                    LineSearchMax = Convert.ToInt16(parts[1]);

                if (line.Contains(names[10]))
                    AlphaMax = Convert.ToDouble(parts[1]);

                if (line.Contains(names[11]))
                    MinCatchBlock = Convert.ToInt16(parts[1]);

                if (line.Contains(names[12]))
                    MinTowingSpeed = Convert.ToDouble(parts[1]);

                if (line.Contains(names[13]))
                    UsePreviousAsPrecalc = Convert.ToBoolean(Convert.ToInt16(parts[1]));

                currentLine++;
            }

        }

        public void PrintInfo()
        {
            Console.WriteLine("\nSolver settings info:");
            Console.WriteLine("{0,-35}{1,-10:D}", "Maximum iteration count", IterMax);
            Console.WriteLine("{0,-35}{1,-10:E3}", "Residual tolerance", ResidualTol);
            Console.WriteLine("{0,-35}{1,-10:E3}", "Maximum residual", ResidualMax);
            Console.WriteLine("{0,-35}{1,-10:E3}", "Displacement tolerance", DisplacementTol);
            Console.WriteLine("{0,-35}{1,-10:E3}", "Added diagonal stiffness", DiagStiffness);
            Console.WriteLine("{0,-35}{1,-10:E3}", "Added stiffness threshold", StiffnessTol);
            Console.WriteLine("{0,-35}{1,-10:E3}", "Added stiffness increase factor", IncreaseStiffnessBy);
            Console.WriteLine("{0,-35}{1,-10:E3}", "Added stiffness reduce factor", ReduceStiffnessBy);
            Console.WriteLine("{0,-35}{1,-10:D}", "Maximum linesearch iteration count", LineSearchMax);
            Console.WriteLine("{0,-35}{1,-10:E3}", "Step limit factor", AlphaMax);
            Console.WriteLine("{0,-35}{1,-10:D}", "Minimum meshes blocked by catch", MinCatchBlock);
            Console.WriteLine("{0,-35}{1,-10:F2}", "Minimum towing speed", MinTowingSpeed);


        }
    }
}
