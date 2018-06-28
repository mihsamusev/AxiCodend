using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AxiCodend
{
    class HexMeshPanelMaterial
    {
        //======================
        // FIELDS
        //======================

        public double MeshSide { get; set; }
        public double KnotSize { get; set; }
        public double TwineEA { get; set; }
        public double KnotEA { get; set; }
        public int MeshOrientation { get; set; }

        #region
        //public double TwineThickness { get; set; }
        //public bool IsDoubleTwine { get; set; }
        //public double InitialOpeningAngle { get; set; }
        //public double MinimumOpeningAngle { get; set; }            
        //public double EI { get; set; }
        //public double OpenningStifness { get; set; }
        //public double KnotContactStifness { get; set; }
        //public double Density { get; set; }
        #endregion
        //======================
        // CONSTRUCTORS
        //======================

        public HexMeshPanelMaterial()
        {
            LoadDefault();
        }

        public HexMeshPanelMaterial(PathsIO path)
        {
            LoadInput(path.input);
        }

        //======================
        // METHODS
        //======================

        private void LoadDefault()
        {
            MeshSide = 0.1;
            KnotSize = 0.001;
            TwineEA = 1000;
            KnotEA = 1000;
            MeshOrientation = 0;
        }

        private void LoadInput(string inpPath)
        {
            string[] names = { "MeshSide", "KnotSize", "TwineEA", "KnotEA", "MeshOrientation" };
            string[] lines = File.ReadAllLines(inpPath);
            string[] parts;
            int currentLine = 0;
            int foundCount = 0;

            foreach (var line in lines)
            {
                if (line.Contains(names[0]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    MeshSide = Convert.ToDouble(parts[1]);
                    foundCount++;
                }

                if (line.Contains(names[1]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    KnotSize = Convert.ToDouble(parts[1]);
                    foundCount++;
                }

                if (line.Contains(names[2]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    TwineEA = Convert.ToDouble(parts[1]);
                    foundCount++;
                }

                if (line.Contains(names[3]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    KnotEA = Convert.ToDouble(parts[1]);
                    foundCount++;
                }

                if (line.Contains(names[4]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    MeshOrientation = Convert.ToInt32(parts[1]);
                    foundCount++;
                }
                currentLine++;
            }

            if(foundCount != 5)
            {
                throw new ArgumentException("Not all fields could be initialized, " +
                                            "because the input file is not in the right format");
            }
        }

        public void PrintInfo()
        {
            Console.WriteLine("\nHexagonal panel material info:");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Mesh side", MeshSide, "[m]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Knot size", KnotSize, "[m]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Twine EA", TwineEA, "[N]");
            Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Knot EA", KnotEA, "[N]");
            Console.WriteLine("{0,-25}{1,-10}", "Orientation", "T"+ MeshOrientation);

            #region
            //Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Twine thickness", TwineThickness, "[m]");
            //Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Initial opening angle", InitialOpeningAngle, "[deg]");
            //if (IsDoubleTwine)
            //{
            //    Console.WriteLine("The mesh is double twine");
            //}
            //else
            //{
            //    Console.WriteLine("The mesh is single twine");
            //}
            //Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Density", Density, "[kg / m3]");

            //Console.WriteLine("{0,-25}{1,-10:F3}{2}", "EI", EI, "[N * m^2]");
            //Console.WriteLine("{0,-25}{1,-10:F3}{2}", "Opening Stiffness", OpenningStifness, "[N * rad]");
            //Console.WriteLine();
            #endregion
        }

        



    }
}
