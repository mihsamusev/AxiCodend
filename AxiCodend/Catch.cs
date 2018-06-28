using System;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AxiCodend
{
    class Catch
    {
        //=========================
        // variables
        //=========================

        public int[] BlockedMeshes;
        public int Count;


        //=========================
        // constructors
        //=========================

        public Catch(int[] BlockedMeshes)
        {
            this.BlockedMeshes = BlockedMeshes;
            Count = BlockedMeshes.Length;
        }

        public Catch(PathsIO path)
        {
            LoadInput(path.input);
        }

        //=========================
        // methods
        //=========================

        private void LoadInput(string inpPath)
        {
            string[] lines = File.ReadAllLines(inpPath);
            int currentLine;
            string[] parts;
            bool catchFound = false;

            // Node List
            currentLine = 0;
            foreach (var line in lines)
            {
                if (line.Contains("Catch") || line.Contains("CATCH"))
                {
                    catchFound = true;
                    break;
                }
                currentLine++;
            }

            if(catchFound)
            {
                currentLine++;
            }
            else
            {
                throw new Exception("Incomplete input file, \'Catch\' or \'CATCH\' tag is not found!");
            }
            
            parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
            Count = Convert.ToInt32(parts[1]);
            if (Count > 0)
            {
                currentLine++;
                BlockedMeshes = new int[Count];
                for (int i = 0; i < Count; i++)
                {
                    parts = lines[currentLine + i].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);   
                    BlockedMeshes[i] = Convert.ToInt32(parts[0]);                   
                }
            }
        }

    }
}
