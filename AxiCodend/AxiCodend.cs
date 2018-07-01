using System;
using System.IO;
using CSparse.Storage;
using CSparse;

namespace AxiCodend
{
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

        /*FEM*/
        public int dof;                    // number of d.o.f.
        public double R;                            // residue
        public double[] X;                          // Shape
        public double[] F;                          // force vectors
        protected CoordinateStorage<double> Jcs;   // Sparse Jacobians
        public CompressedColumnStorage<double> J;

        /*others*/
        private const double rhoWater = 1025;
        private const double catchCd = 1.4;
        //====================
        // CLASS CONSTRUCTOR
        //====================

        public AxiCodend()
        {

        }

        public AxiCodend(int nx, int nr, double r0, HexMeshPanelMaterial Material)
        {
            /*main inputs*/
            this.nx = nx;
            this.nr = nr;
            this.r0 = r0;
            this.Material = Material;

            nc = nx;    // to start with

            l0 = Material.KnotSize;
            m0 = Material.MeshSide / 2;
            kl = Material.KnotEA;
            km = Material.TwineEA;

            InitP();
        }

        public AxiCodend(PathsIO path)
        {
            LoadInput(path);
            InitP();
        }

        //====================
        // CLASS METHODS
        //====================

        /*private methods*/

        private void LoadInput(PathsIO path)
        {

        }

        private void InitP()
        {
            towSpeed = 1;
            P = 0.5 * catchCd * rhoWater * Math.Pow(towSpeed, 2);
        }

        protected virtual void InitDOF()
        {

        }

        protected virtual void InitBlockedMeshes()
        {

        }

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

        public virtual void ClearState()
        {

        }

        public virtual void ApplyCatch(int newCatch)
        {

        }

        public virtual void ApplyTowing(double newSpeed)
        {

        }

        public virtual void UpdateTotalJacobian(double kDiag, bool IncludeBC)
        {

        }

        public virtual void UpdatePosition(double[] h, double lambda)
        {

        }

        public virtual void UpdateTotalForces(bool IncludeBC)
        {

        }

        public virtual void UpdateResidual()
        {

        }

        public virtual void SetInitialShape()
        {

        }

        public virtual void SetInitialShape(string path)
        {

        }

        public virtual void SetInitialShapeSmooth()
        {

        }



        public virtual double Length()
        {
            return 0;
        }

        public virtual double MaxRadius()
        {
            return 0;
        }

        public virtual double EntranceDrag()
        {
            return 0;
        }

        public virtual double CatchDrag()
        {
            return 0;
        }

        public virtual double CatchThickness()
        {
            return 0;
        }

        public virtual double CatchVolume()
        {
            return 0;
        }

        public virtual double CatchSurface()
        {
            return 0;
        }
    }
}
