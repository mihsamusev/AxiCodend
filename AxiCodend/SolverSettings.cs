namespace AxiCodend
{
    struct SolverSettings
    {
            public int IterMax { get; set; } = 5000;
            public double ResidualTol { get; set; } = 1e-3;                  
            public double DisplacementTol { get; set; } = 1e-4;                 
            public double ResidualMax { get; set; } = 10e20;
            public double StiffnessTol { get; set; } = 1.0;
            public double DiagStiffness { get; set; } = 1.0;
            public double ReduceStiffnessBy { get; set; } = 0.1;
            public double IncreaseStiffnessBy { get; set; } = 2.0;
            public bool ShowLineSearchSteps { get; set; } = false;
            public int LineSearchMax { get; set; } = 6;
            public double AlphaMax { get; set; } = 0.2 * 0.16; 
            public int MinCatchBlock { get; set; } = 5;
            public double MinTowingSpeed { get; set; } = 0.1;
            public bool UsePreviousAsPrecalc { get; set; } = false;

        public override string ToString()
        {
            return String.Format("\nSolver settings info:\n") +
                String.Format("{0,-35}{1,-10:D}\n", "Maximum iteration count", IterMax) +
                String.Format("{0,-35}{1,-10:E3}\n", "Residual tolerance", ResidualTol) +
                String.Format("{0,-35}{1,-10:E3}\n", "Maximum residual", ResidualMax) +
                String.Format("{0,-35}{1,-10:E3}\n", "Displacement tolerance", DisplacementTol) +
                String.Format("{0,-35}{1,-10:E3}\n", "Added diagonal stiffness", DiagStiffness) +
                String.Format("{0,-35}{1,-10:E3}\n", "Added stiffness threshold", StiffnessTol) +
                String.Format("{0,-35}{1,-10:E3}\n", "Added stiffness increase factor", IncreaseStiffnessBy) +
                String.Format("{0,-35}{1,-10:E3}\n", "Added stiffness reduce factor", ReduceStiffnessBy) +
                String.Format("{0,-35}{1,-10:D}\n", "Maximum linesearch iteration count", LineSearchMax) +
                String.Format("{0,-35}{1,-10:E3}\n", "Step limit factor", AlphaMax) +
                String.Format("{0,-35}{1,-10:D}\n", "Minimum meshes blocked by catch", MinCatchBlock) +
                String.Format("{0,-35}{1,-10:F2}\n", "Minimum towing speed", MinTowingSpeed);
        }
    }
}
