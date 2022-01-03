using System;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using CSparse.Storage;
using CSparse;
using CSparse.Double.Factorization;

namespace AxiCodend
{
    class Program
    {
        static void Main(string[] args)
        {
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
                codend = new AxiModelT0(path);   // initialize T0 model
            }
            else if(hexMeshPanel.MeshOrientation == 90)
            {
                codend = new AxiModelT90(path);   // initialize T90 model
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
        }


    }
}