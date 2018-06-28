using System;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AxiCodend
{
    class Towing
    {
        //=========================
        // variables
        //=========================

        public double Speed;

        //=========================
        // constructor
        //=========================

        public Towing(double Speed)
        {
            this.Speed = Speed;
        }

        public Towing(PathsIO path)
        {
            LoadInput(path.input);
        }

        //=========================
        // methods
        //=========================

        private void LoadInput(string inpPath)
        {
            string[] names = { "TowingSpeed" };
            string[] lines = File.ReadAllLines(inpPath);
            string[] parts;
            int currentLine = 0;

            foreach (var line in lines)
            {
                if (line.Contains(names[0]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    Speed = Convert.ToDouble(parts[1]);
                    break;
                }
                currentLine++;
            }
        }

    }
}
