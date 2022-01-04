

using CommandLine;
using YamlDotNet.Serialization.NamingConventions;

namespace AxiCodend
{
    struct CodendGeometry {
        public int MeshesAlong {get; set;} = 100;                 
        public int MeshesAround {get; set;} = 100;                 
        public double EntranceRadius {get; set;} = 0.400000;

        public override string ToString()
        {
            return "Codend geometry:\n" +  
                String.Format("{0,-35}{1,-10:D}\n", "Meshes along", MeshesAlong) + 
                String.Format("{0,-35}{1,-10:D}\n", "Meshes around", MeshesAround) + 
                String.Format("{0,-35}{1,-10:F2}{2}\n", "Entrance radius", EntranceRadius,"[m]");
        }
    }

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
            Console.WriteLine(hexMeshMaterial.ToString());

            var solverSettings = cfg.solver;
            Console.WriteLine(solverSettings.ToString());
            
            var codendGeometry = cfg.geometry;
            Console.WriteLine(codendGeometry.ToString());
            
            AxiCodend codend = new AxiCodend();
            
            if (hexMeshMaterial.MeshOrientation == MeshOrientation.T0)
            {
                codend = new AxiModelT0(codendGeometry, hexMeshMaterial);   // initialize T0 model
            }
            else
            {
                throw new IOException("Incorrect mesh orientation. Value should be 0 or 90.");
            }

            var TowingSimulation = new Simulation(codend, cfg.catches, cfg.towing_speed)
            {
                SolverSettings = solverSettings,
                Paths = cfg.paths
            };

            TowingSimulation.Simulate();
        }


    }
}