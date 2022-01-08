using CommandLine;
using YamlDotNet.Serialization.NamingConventions;

namespace AxiCodend
{

    enum OutputFormats {
        txt = 0,
        json = 1
    }

    struct Output {
        public string Filename {get; init;} = "results";
        public OutputFormats Format {get; init;} = OutputFormats.json;

        public bool Valid() {
            return Directory.Exists(Path.GetDirectoryName(Filename));
        }

        public string GetPath() {
            return Path.Combine(
                Directory.GetCurrentDirectory(), 
                Path.ChangeExtension(this.Filename, this.Format.ToString()));
        }
    }

    public class InputArguments
    {
        [Option('j', "job", Required = true, HelpText = "YAML job configuration file")]
        public string? jobFile {get; set;}
    }

    struct JobConfig {
        public Output output {get; set;}
        public HexMeshPanelMaterial material {get; set;}
        public CodendGeometry geometry {get; set;}
        public int[] catches {get; set;} = {};
        public double towing_speed {get; set;} = 0;
        public SolverSettings solver {get; set;}

        public JobConfig() {
            output = new Output();
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
            var codendMaterial = cfg.material;
            var solverSettings = cfg.solver;
            var codendGeometry = cfg.geometry;
            

            AxiCodend codend = new AxiCodend();
            if (codendMaterial.MeshOrientation == MeshOrientation.T0) {
                codend = new AxiModelT0(codendGeometry, codendMaterial);   // initialize T0 model
            } else if (codendMaterial.MeshOrientation == MeshOrientation.T90) {
                codend = new AxiModelT90(codendGeometry, codendMaterial);
            } else {
                throw new ArgumentException(
                    "Incorrect mesh orientation. Value should be 0 or 90.");
            }

            var output = cfg.output;
            ICodendSaver resultSaver;
            if (output.Format == OutputFormats.json) {
                resultSaver = new JSONResultSaver(output.GetPath());
            } else if (output.Format == OutputFormats.txt) {
                resultSaver = new CSVResultSaver(output.GetPath());
            } else {
                throw new ArgumentException(
                    "Incorrect file format. Value should \"json\" or \"txt\".");  
            }

            var TowingSimulation = new Simulation(
                codend, cfg.catches, cfg.towing_speed, resultSaver) {
                SolverSettings = solverSettings,
            };

            TowingSimulation.Simulate();
        }


    }
}