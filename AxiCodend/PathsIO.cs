using System;
using System.Reflection;
using System.IO;

namespace AxiCodend
{
    class PathsIO
    {
        //=========================
        // fields
        //=========================

        public string input;
        public string outputShapes;
        public string outputResults;

        //=========================
        // constructors
        //=========================

        public PathsIO()
        {    
            string dir = "/home/msa/Documents/se_path/cs/AxiCodend/AxiCodend/bin/Debug/net6.0/";
            input = Path.Combine(dir, "input.txt");
            outputShapes = Path.Combine(dir, "shape.txt");
            outputResults = Path.Combine(dir, "results.txt");
        }

        //=========================
        // methods
        //=========================


        public void PrintInfo()
        {
            Console.WriteLine("\nPATHS:");
            Console.WriteLine("Input to the current 3d model:");
            Console.WriteLine("\t" + input);
            Console.WriteLine("Output shape from axis-symmetric model:");
            Console.WriteLine("\t" + outputShapes);
            Console.WriteLine("Output result parameters from the current 3d model:");
            Console.WriteLine("\t" + outputResults);
            Console.WriteLine();
        }

    }
}
