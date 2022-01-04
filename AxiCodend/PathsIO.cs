using System;
using System.Reflection;
using System.IO;

namespace AxiCodend
{
    struct PathsIO
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


    struct OutputPaths
    {
        public string OutputShapes {get; set;}
        public string OutputResults {get; set;}

        public bool Valid() {
            var isValidShapesDir = Directory.Exists(Path.GetDirectoryName(OutputShapes));
            var isValidResultsDir = Directory.Exists(Path.GetDirectoryName(OutputResults));
            return isValidShapesDir && isValidResultsDir;
        }

        public OutputPaths() {
            OutputShapes = Path.Combine(Directory.GetCurrentDirectory(), "shapes.txt");
            OutputResults = Path.Combine(Directory.GetCurrentDirectory(), "results.txt");
        }

        public override string ToString() {
            return "\nResult save path:\n" + 
                String.Format("Simulated shapes are saved to {0}", OutputShapes) + 
                String.Format("Simulated results are saved to {0}", OutputResults);
        } 
    }
}
