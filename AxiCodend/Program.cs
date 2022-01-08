

using CommandLine;
using YamlDotNet.Serialization.NamingConventions;

namespace AxiCodend
{


    public class InputArguments
    {
        [Option('j', "job", Required = true, HelpText = "YAML job configuration file")]
        public string? jobFile {get; set;}
    }

    struct JobConfig
    {
        public OutputPaths paths {get; set;}
        public HexMeshPanelMaterial material {get; set;}
        public CodendGeometry geometry {get; set;}
        public int[] catches {get; set;} = {};
        public double towing_speed {get; set;} = 0;
        public SolverSettings solver {get; set;}

        public JobConfig() {
            paths = new OutputPaths();
            material = new HexMeshPanelMaterial();
            geometry = new CodendGeometry();
            solver = new SolverSettings();
        }
    }

    class Program
    {
        static T ParseYaml<T>(string configPath)
        {
        var deserializer = new YamlDotNet.Serialization.DeserializerBuilder()
            .WithNamingConvention(UnderscoredNamingConvention.Instance)
            .Build();

        return deserializer.Deserialize<T>(
            File.ReadAllText(configPath));
        }

        static void Main(string[] args) {
            var results = Parser.Default.ParseArguments<InputArguments>(args)
                .WithParsed(args => Run(args));
        }


        static void Run(InputArguments args)
        {
            var cfg = ParseYaml<JobConfig>(args.jobFile);
            var hexMeshMaterial = cfg.material;
            var solverSettings = cfg.solver;
            var codendGeometry = cfg.geometry;

            
            AxiCodend codend = new AxiCodend();
            if (hexMeshMaterial.MeshOrientation == MeshOrientation.T0) {
                codend = new AxiModelT0(codendGeometry, hexMeshMaterial);   // initialize T0 model
            } else if (hexMeshMaterial.MeshOrientation == MeshOrientation.T90) {
                codend = new AxiModelT90(codendGeometry, hexMeshMaterial);
            } else {
                throw new IOException("Incorrect mesh orientation. Value should be 0 or 90.");
            }

            var rs = new CSVResultSaver(cfg.paths.OutputShapes);

            // var csvSaver = new CSVResultSaverBuilder
            //    .WithHeader(hexMeshMaterial, codendGeometry);
            // var jsonSaver = new JSONResultSaverBuilder
            //    .WithHeader(hexMeshMaterial, codendGeometry);

            var TowingSimulation = new Simulation(codend, cfg.catches, cfg.towing_speed, rs)
            {
                SolverSettings = solverSettings,
                Paths = cfg.paths
            };

            TowingSimulation.Simulate();

            Console.Write(codend.GetMetrics().ToString());
        }


    }
}