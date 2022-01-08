

namespace AxiCodend
{
    public enum MeshOrientation
    {
        T0 = 0,
        T90 = 90
    }

    public struct HexMeshPanelMaterial
    {
        //======================
        // FIELDS
        //======================

        public double MeshSide { get; set; } = 0.1;
        public double KnotSize { get; set; } = 0.001;
        public double TwineStiffness { get; set; } = 1000;
        public double KnotStiffness { get; set; } = 1000;
        public MeshOrientation MeshOrientation { get; set; } = 0;


        public HexMeshPanelMaterial(){}

        public override string ToString()
        {
            return String.Format("Hexagonal panel material info:\n") + 
                String.Format("{0,-35}{1,-10:F3}{2}\n", "Mesh side", MeshSide, "[m]") +
                String.Format("{0,-35}{1,-10:F3}{2}\n", "Knot size", KnotSize, "[m]") +
                String.Format("{0,-35}{1,-10:F3}{2}\n", "Twine tensile stiffness [EA]", TwineStiffness, "[N]") +
                String.Format("{0,-35}{1,-10:F3}{2}\n", "Knot tensile stiffness [EA]", KnotStiffness, "[N]") +
                String.Format("{0,-35}{1,-10}\n", "Orientation", MeshOrientation);
        }
    }
}
