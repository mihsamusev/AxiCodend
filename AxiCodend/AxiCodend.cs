using System;
using System.IO;
using CSparse.Storage;
using CSparse;

namespace AxiCodend
{
    public struct CodendGeometry {
        public int MeshesAlong {get; set;} = 100;                 
        public int MeshesAround {get; set;} = 100;                 
        public double EntranceRadius {get; set;} = 0.400000;

        public override string ToString()
        {
            return "Codend geometry:\n" +  
                String.Format("{0,-35}{1,-10:D}\n", "Meshes along", MeshesAlong) + 
                String.Format("{0,-35}{1,-10:D}\n", "Meshes around", MeshesAround) + 
                String.Format("{0,-35}{1,-10:F2}{2}\n", "Entrance radius", EntranceRadius,"[m]");
        }
    }


    public struct CodendMetrics {
        public double Length {get; init;}
        public double MaxRadius {get; init;}
        public double CatchThickness {get; init;}
        public double CatchSurface {get; init;}
        public double CatchVolume {get; init;}
        public double CatchDrag {get; init;}
        public double EntranceDrag {get; init;}

        public override string ToString()
        {
            return String.Format("Metrics:\n") + 
                String.Format("Total codend length:             {0,10:N3} [m]\n", Length) +
                String.Format("Maximum codend radius:           {0,10:N3} [m]\n", MaxRadius) +
                String.Format("Catch extent:                    {0,10:N3} [m]\n", CatchThickness) +
                String.Format("Surface in contact with catch:   {0,10:N3} [m^2]\n", CatchSurface) +
                String.Format("Catch volume:                    {0,10:N3} [m^3]\n", CatchVolume) + 
                String.Format("Resultant entrance reaction:     {0,10:N0} [N]\n", EntranceDrag) +
                String.Format("Total catch drag force:          {0,10:N0} [N]\n", CatchDrag);
        }
    }


    class AxiCodend
    {
        //====================
        // CLASS VARIABLES
        //====================

        /*codend parameters*/
        protected double l0;           // length of un-stretched knot       
        protected double m0;           // length of un-stretched twine
        protected double kl;           // stiffness of a knot
        protected double km;           // stiffness of a twine

        protected int nx;              // amount of meshes along the codend length
        protected int nr;              // amount of meshes in the codend`s circumference
        protected int nc;              // moundt of meshes along the codend length contatining catch 
        protected double r0;           // entrance radius
        protected int ncp;             // node number indicating start of application of the catch pressure

        /*operation*/
        public double towSpeed;        // towing speed
        protected double P;            // catch pressure

        /*material*/
        public HexMeshPanelMaterial Material;
        public CodendGeometry Geometry;

        /*FEM*/
        public int dof;                    // number of d.o.f.
        public double R;                            // residue
        public double[] X;                          // Shape
        public double[] F;                          // force vectors
        protected CoordinateStorage<double> Jcs;   // Sparse Jacobians
        public CompressedColumnStorage<double> J;

        /*others*/
        protected const double rhoWater = 1025;
        protected const double catchCd = 1.4;

        //====================
        // CLASS CONSTRUCTOR
        //====================

        public AxiCodend()
        {

        }

        public AxiCodend(CodendGeometry geometry, HexMeshPanelMaterial material)
        {
            /*main inputs*/
            this.Geometry = geometry;
            this.Material = material;

            this.nx = geometry.MeshesAlong;
            this.nr = geometry.MeshesAround;
            this.r0 = geometry.EntranceRadius;
            
            nc = nx;    // to start with entire codend is blocked

            l0 = material.KnotSize;
            m0 = material.MeshSide / 2;
            kl = material.KnotStiffness;
            km = material.TwineStiffness;

            towSpeed = 1;
            P = 0.5 * catchCd * rhoWater * Math.Pow(towSpeed, 2);

            SetDOF();
            SetBlockedMeshes();
            ClearState();
            SetAngles();
        }

        //====================
        // CLASS METHODS
        //====================

        /*private methods*/

        private void SetPressure()
        {
            towSpeed = 1;
            P = 0.5 * catchCd * rhoWater * Math.Pow(towSpeed, 2);
        }

        protected virtual void SetDOF() {}

        protected virtual void SetBlockedMeshes() {}

        protected virtual void SetAngles() {}

        /* public methods*/

        public int GetMeshesBlocked()
        {
            return nc;
        }

        public int GetMeshesAlong()
        {
            return nx;
        }

        public int GetMeshesAround()
        {
            return nr;
        }

        public virtual void ClearState() {}

        public virtual void ApplyCatch(int newCatch) {}

        public virtual void ApplyTowing(double newSpeed) {}

        public virtual void UpdateTotalJacobian(double kDiag, bool IncludeBC) {}

        public virtual void UpdatePosition(double[] h, double lambda) {}

        public virtual void UpdateTotalForces(bool IncludeBC) {}

        public virtual void UpdateResidual() {}

        public virtual void SetInitialShape() {}

        public virtual void SetInitialShape(string path) {}

        public virtual void SetInitialShapeSmooth() {}

        public virtual double Length() {return 0;}

        public virtual double MaxRadius()  {return 0;}

        public virtual double EntranceDrag() {return 0;}

        public virtual double CatchDrag()  {return 0;}

        public virtual double CatchThickness() {return 0;}

        public virtual double CatchVolume() {return 0;}

        public virtual double CatchSurface() {return 0;}

        public CodendMetrics GetMetrics()
        {
            return new CodendMetrics() {
                Length = Length(),
                MaxRadius = MaxRadius(),
                CatchThickness = CatchThickness(),
                CatchSurface = CatchSurface(),
                CatchVolume = CatchVolume(),
                EntranceDrag = EntranceDrag(),
                CatchDrag = CatchDrag()
            };
        }
    }
}
