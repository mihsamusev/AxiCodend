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
        private double l0;           // length of un-stretched knot       
        private double m0;           // length of un-stretched twine
        private double kl;           // stiffness of a knot
        private double km;           // stiffness of a twine
   
        private int nx;              // amount of meshes along the codend length
        private int nr;              // amount of meshes in the codend`s circumference
        private int nc;              // moundt of meshes along the codend length contatining catch 
        private double r0;           // entrance radius

        /*operation*/
        public double towSpeed;     // towing speed
        private double P;           // catch pressure

        /*material*/
        public HexMeshPanelMaterial Material;

        /*FEM*/
        public readonly int dof;                 // number of d.o.f.
        public double R;                         // residue
        public double[] X;                       // Shape
        public double[] F;                       // force vectors
        private CoordinateStorage<double> Jcs;   // Sparse Jacobians
        public CompressedColumnStorage<double> J;

        /*others*/
        private const double rhoWater = 1025;
        private const double catchCd = 1.4;
        //====================
        // CLASS CONSTRUCTOR
        //====================

        public AxiCodend(int nx, int nr, double r0, HexMeshPanelMaterial Material)
        {
            /*main inputs*/
            this.nx = nx;
            this.nr = nr;
            this.r0 = r0;
            this.Material = Material;

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

        private void LoadInput(PathsIO path)
        {

        }

        private void InitP()
        {
            towSpeed = 1;
            P = 0.5 * catchCd * rhoWater * Math.Pow(towSpeed, 2);
        }

    }
}
