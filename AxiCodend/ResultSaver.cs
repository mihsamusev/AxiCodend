using System.IO;
using System.Text;

namespace AxiCodend
{
    interface ICodendSaver {
        string Path {get; set;}
        void append_run(int blockedMeshes, double towingSpeed, CodendMetrics metrics, double[] dofShape);
        void header(HexMeshPanelMaterial material, CodendGeometry geometry);
        void save(double[][] shapes);
        // void append(double[] shape);
    }


    public class CSVResultSaver : ICodendSaver
    {
        public string Path {get; set;}
        public string Separator {get; set;}
        private int HeaderSize = 27;

        public CSVResultSaver(string path, string separator = "\t") {
            this.Path = path;
            this.Separator = separator;
        }

        public void header(HexMeshPanelMaterial material, CodendGeometry geometry) {
            using (StreamWriter ResultFile = new StreamWriter(Path)) {
                ResultFile.WriteLine("HEADER:");
                ResultFile.WriteLine(material.ToString());
                ResultFile.WriteLine(geometry.ToString());
                ResultFile.WriteLine("COLUMN FORMAT:");
                ResultFile.WriteLine("Applied Catch");
                ResultFile.WriteLine("Towing Speed");
                ResultFile.WriteLine("Length");
                ResultFile.WriteLine("MaxRadius");
                ResultFile.WriteLine("CatchThickness");
                ResultFile.WriteLine("CatchSurface");
                ResultFile.WriteLine("CatchVolume");
                ResultFile.WriteLine("CatchDrag");
                ResultFile.WriteLine("EntranceDrag");
                ResultFile.WriteLine("DofShape\n");
                ResultFile.WriteLine("RUNS:\n");
            }
        }

        public void append_run(int blockedMeshes, double towingSpeed, CodendMetrics metrics, double[] dofShape) {
            string[] lines = File.ReadAllLines(Path);
            string output = "";
            if (lines.Length > HeaderSize) {
                var sb = new StringBuilder(
                    string.Join("\n", lines, 0, HeaderSize -1) + "\n",
                    2 * dofShape.Length
                );

                sb.AppendLine(lines[HeaderSize - 1] + Separator +  String.Format("{0:N2}", blockedMeshes));
                sb.AppendLine(lines[HeaderSize    ] + Separator +  String.Format("{0:N2}", towingSpeed));
                sb.AppendLine(lines[HeaderSize + 1] + Separator +  String.Format("{0:N2}", metrics.Length));
                sb.AppendLine(lines[HeaderSize + 2] + Separator +  String.Format("{0:N2}", metrics.MaxRadius));
                sb.AppendLine(lines[HeaderSize + 3] + Separator +  String.Format("{0:N2}", metrics.CatchThickness));
                sb.AppendLine(lines[HeaderSize + 4] + Separator +  String.Format("{0:N2}", metrics.CatchVolume));
                sb.AppendLine(lines[HeaderSize + 5] + Separator +  String.Format("{0:N2}", metrics.CatchDrag));
                sb.AppendLine(lines[HeaderSize + 6] + Separator +  String.Format("{0:N2}", metrics.EntranceDrag));

                // append dof shape
                for (int i = 0; i < dofShape.Length; i++) {
                    sb.AppendLine(lines[HeaderSize + 7] + Separator + String.Format("{0:E5}", dofShape[i]));
                }
                output = sb.ToString();
            } else {
                var sb = new StringBuilder(
                    string.Join("\n", lines),
                    2 * dofShape.Length
                );

                // append metrics
                sb.AppendLine(String.Format("{0:N2}", blockedMeshes));
                sb.AppendLine(String.Format("{0:N2}", towingSpeed));
                sb.AppendLine(String.Format("{0:N2}", metrics.Length));
                sb.AppendLine(String.Format("{0:N2}", metrics.MaxRadius));
                sb.AppendLine(String.Format("{0:N2}", metrics.CatchThickness));
                sb.AppendLine(String.Format("{0:N2}", metrics.CatchVolume));
                sb.AppendLine(String.Format("{0:N2}", metrics.CatchDrag));
                sb.AppendLine(String.Format("{0:N2}", metrics.EntranceDrag));

                // append dof shape
                foreach(var dof in dofShape) {
                    sb.AppendLine(String.Format("{0:E5}", dof));
                }
                output = sb.ToString();
            }
            File.WriteAllText(Path, output);
        }

        public void save(double[][] shapes) {
            using (StreamWriter ResultFile = File.AppendText(Path))
            {
                for (int row = 0; row < shapes[0].Length; row++)
                {
                    for (int col = 0; col < shapes.Length - 1; col++)
                    {
                        ResultFile.Write("{0:E5}{1}", shapes[col][row], this.Separator);
                    }
                    ResultFile.Write("{0:E5}\n", shapes[shapes.Length - 1][row]);
                }
            }
        }
        // end CSVShapeSaver
    }

    // public class JSONResultSaver : ICodendSaver {
    //   string Path
    //}

}
