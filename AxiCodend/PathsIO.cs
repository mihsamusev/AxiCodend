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
            LoadDefault();
            LoadFromInput(input);
        }

        //=========================
        // methods
        //=========================

        public string AssemblyPath()
        {
            string fullPath = Path.GetDirectoryName(Assembly.GetExecutingAssembly().GetName().CodeBase);
            return fullPath.Replace("file:\\", "");
        }

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

        private void LoadDefault()
        {
            input = AssemblyPath() + "\\input.txt";
            outputShapes = AssemblyPath() + "\\shape.txt";
            outputResults = AssemblyPath() + "\\results.txt";
        }

        private void LoadFromInput(string input)
        {
            string[] names = {"OutputShapes", "OutputResults" };
            string[] lines = File.ReadAllLines(input);
            string[] parts;
            int currentLine = 0;

            foreach (var line in lines)
            {
                parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);

                if (line.Contains(names[0]))
                    outputShapes = parts[1];

                if (line.Contains(names[1]))
                    outputResults = parts[1];

                currentLine++;
            }
        }
}
}
