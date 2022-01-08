# AxiCodend
A cod-end is the rearmost part of a trawl fishing gear that
collects the catch during the towing. The shape of cod-ends
is of importance as it determines mesh opening and consequently
influences the selectivity of fish from the cod-end. Poor selectivity results in negative environmental impact due to collection of large quantities of juvenile fish that is usually discarded as by-catch. `AxiCodend` is a CLI tool for simulation of trawl cod-end shapes. It is developed as a object-oriented C# version of axis-symmetric numerical model developed by D.Priour [[Paper](https://www.sciencedirect.com/science/article/abs/pii/S0029801814003709?via%3Dihub)]. This implementation was created as a part of the master thesis in collaboration with SINTEF Ocean: _Implementation and comparison of two numerical models of trawl cod-end_ [[Link](https://projekter.aau.dk/projekter/en/studentthesis/implementation-and-comparison-of-two-numerical-models-for-trawl-codends(7c4900a9-f83e-4f61-818b-2c271252cab1).html)].

![](/docs/codend_example.png)

_Note: this visualization is not included_.

## How to build
- Install [.NET 6.0 Core SDK](https://dotnet.microsoft.com/en-us/download)
- Install [git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)

Using Terminal or PowerShell:
```sh
git clone https://github.com/mihsamusev/AxiCodend \
cd AxiCodend \
dotnet build --configuration Release
```

## Getting started
A simulation job can be run using the `dotnet` command

```sh
dotnet run --project AxiCodend --job example_job.yaml # from the solution folder
```

The job configuration file consists of settings regarding the cod-end geometry, materials, loads, simulation output paths and solver settings.

```yaml
material: 
  mesh_side: 0.120000           # twine length in [m]
  mesh_orientation: 0           # T0 or T90 (see paper)
  knot_size: 0.012000           # thickness of the knot [m]
  twine_stiffness: 1000.000000  # axial stiffness EA in [Nm2]
  knot_stiffness: 1000.000000   # axial stiffness EA in [Nm2]

geometry: 
  meshes_along: 100                 
  meshes_around: 100                 
  entrance_radius: 0.400000 # [m]

catches: [50, 75]           # in number of blocked meshes from the end
towing_speed: 1.500000      # speed [m/s]

output: 
  filename: myresults
  format: json                  # [json, txt]
```

Optionally, one can configure many the iterative solver parameters. The method is based on Newton-Raphson method with a few additional heuristics for descent stability. Based on our testing the solver parameters given below are sufficient for 99% of the cases. However if in need of more fine grained control the simulation algorithm and solver parameters are described in the thesis.

```yaml
solver:
  iter_max: 3000
  residual_tol: 1e-3         # when to finishe the iteration based on norm force < 0.001 N
  displacement_tol: 1e-4     # when to finishe the iteration based on norm dispacement < 0.1 mm
  residual_max: 10e20        # Maximum residual after which the scheme is considered divergent
  stiffness_tol: 1           #
  diag_stiffness: 1          # additional stiffness for stable iterations

  reduce_stiffness_by: 0.1   # if system diverged due to being too stiff what to do on restart calculation
  increase_stiffness_by: 2   # if sysyem diverged due to being too soft 

  show_line_search_steps: false # show intermedite (line search) Newton method solver steps
  line_search_max: 6            # maximum amount of inexact line search iterations
  alpha_max: 0.032              # alpha for Armijo step limit globalization method (20% of twine length)

  min_catch_block: 5            
  min_towing_speed: 0.1            
  use_previous_as_precalc: false # warm start for the series of simulations
```


### Typical results
If the simulation is successful `.json` result contains a header dublicating the material and the geometry from the job config. The `runs` key contains a list of simulation results, one for each `meshes_blocked` input. The `metrics` contains the general description of the resulting shape like length and volume. the `dof_shape` contains the list of coordinates for each degree of freedom in the resulting codend shape. Depending on the codend mesh orientation the coordinates are polar `[x_0, r_0, ..., x_n, r_n]` for `T0` or cartesian `[x_0, y_0, z_0, ..., x_n, y_n, z_n]` for `T90`. All units are SI.

```json
{
  "material": {
    "mesh_side": 0.12,
    "knot_size": 0.012,
    "twine_stiffness": 1000.0,
    "knot_stiffness": 1000.0,
    "mesh_orientation": 0
  },
  "geometry": {
    "meshes_along": 100,
    "meshes_around": 100,
    "entrance_radius": 0.4
  },
  "runs": [
    {
      "meshes_blocked": 50,
      "towing_speed": 1.5,
      "metrics": {
        "length": 11.677690629036436,
        "max_radius": 1.8757983780600964,
        "catch_thickness": 4.544390246871821,
        "catch_surface": 53.75000030810305,
        "catch_volume": 41.29709565785832,
        "catch_drag": 10069.469117747492,
        "entrance_drag": 10069.469135350682
      },
      "dof_shape": [
        0.0,
        0.4,
        0.013094020279622212,
        0.3981795816161083,
        ...
        11.677690629036436,
        0.0
      ]
    },
    {
      "meshes_blocked": 75,
      "towing_speed": 1.5,
      "metrics": {
        "length": 10.44678264074785,
        "max_radius": 1.9306513710826414,
        "catch_thickness": 7.081715492570243,
        "catch_surface": 84.84607570796004,
        "catch_volume": 71.40537314715513,
        "catch_drag": 10193.596462004796,
        "entrance_drag": 10193.596482952864
      },
      "dof_shape": [
        0.0,
        0.4,
        0.013197741564606073,
        0.40086159222349127,
        0.0748995701004251,
        0.4050895789350485,
        ...
        10.44678264074785,
        0.0
      ]
    }
  ]
}
```
## Visualization
The `preview_results.py` is a small visualization `python` script that depends on `matplotlib` and `numpy` libraries. Run as follows.

```sh
python preview_shapes -i myresults.json
```


## Dependencies
- [Accord.Math](https://www.nuget.org/packages/Accord.Math/) - dense matrix core FEM computations (will be removed asap)
- [CSparse](https://www.nuget.org/packages/CSparse/) - sparse matrix core FEM computations
- [YamlDotNet](https://www.nuget.org/packages/YamlDotNet/) - YAML job configuration parsing
- [CommandLineParser](https://www.nuget.org/packages/CommandLineParser) - CLI argments parsing
- [Newtonsoft.Json](https://www.nuget.org/packages/Newtonsoft.Json/) - JSON Output