using System.Text;
using Newtonsoft.Json;
using Newtonsoft.Json.Serialization;

namespace AxiCodend
{
    interface ICodendSaver {
        string Path {get; init;}
        void AddRun(int blockedMeshes, double towingSpeed, CodendMetrics metrics, double[] dofShape);
        void AddHeader(HexMeshPanelMaterial material, CodendGeometry geometry);
        void Finish();
    }

    public class CSVResultSaver : ICodendSaver
    {
        public string Path {get; init;}
        public string Separator {get; set;}
        private int HeaderSize = 27;

        public CSVResultSaver(string path, string separator = "\t") {
            this.Path = path;
            this.Separator = separator;
        }

        public void AddHeader(HexMeshPanelMaterial material, CodendGeometry geometry) {
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

        public void AddRun(int blockedMeshes, double towingSpeed, CodendMetrics metrics, double[] dofShape) {
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

        public void Finish() {

        }
        // end CSVShapeSaver
    }


    public struct JSONResultTemplate {
        public HexMeshPanelMaterial material;
        public CodendGeometry geometry;
        public List<JSONResultRun> runs = new List<JSONResultRun>();
    }

    public struct JSONResultRun {
        public int meshes_blocked {get; init;}
        public double towing_speed {get; init;}
        public CodendMetrics metrics {get; init;}
        public double[] dof_shape {get; init;}
    }

    public class JSONResultSaver : ICodendSaver {
        public string Path {get; init;}
        private JSONResultTemplate template;

        public JSONResultSaver(string Path) {
            this.Path = Path;
            this.template = new JSONResultTemplate();
        }

        public void AddRun(int blockedMeshes, double towingSpeed, CodendMetrics metrics, double[] dofShape) {
            this.template.runs.Add(
                new JSONResultRun {
                    meshes_blocked = blockedMeshes,
                    towing_speed = towingSpeed,
                    metrics = metrics,
                    dof_shape = dofShape
                }
            );
        }
        
        public void AddHeader(HexMeshPanelMaterial material, CodendGeometry geometry) {
            this.template.material = material;
            this.template.geometry = geometry;
        }

        public void Finish() {
            // based on
            // https://www.newtonsoft.com/json/help/html/NamingStrategySnakeCase.htm
            //
            DefaultContractResolver contractResolver = new DefaultContractResolver
            {
                NamingStrategy = new SnakeCaseNamingStrategy()
            };

            string json = JsonConvert.SerializeObject(
                this.template, new JsonSerializerSettings {
                    ContractResolver = contractResolver,
                    Formatting = Formatting.Indented
            });
            using (StreamWriter ResultFile = new StreamWriter(Path)) {
                ResultFile.Write(json);
            }
        }
    }

}
