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
dotnet build --release
```

## Getting started
A simulation job can be run using the `dotnet` command

```sh
dotnet run -p AxiCodend --job example_job.yaml # from the solution folder
```

The job configuration file consists of settings regarding the cod-end geometry, materials, loads, simulation output paths and solver settings.

```yaml
paths: 
  output_shapes: demoshapes.txt
  output_results: demoresults.txt        

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
*Note: Based on our testing the given solver parameters are sufficient for 99% of the cases. However if in need of more fine grained control the simulation algorithm and solver parameters are described in the thesis.*

### Typical results
If the simulation is successful the results are saved in 2 `.txt` files. The `results.txt` containts global simulation results such as total drag force, achived length and max cod-end radius. 
```
Length              Max radius          Catch thickness     Catch surface       Catch volume        Total reaction      Total force         
1.16777E+001        1.87580E+000        4.54439E+000        5.37500E+001        4.12971E+001        1.00695E+004        1.00695E+004        
1.04468E+001        1.93065E+000        7.08172E+000        8.48461E+001        7.14054E+001        1.01936E+004        1.01936E+004   
```

The `shapes.txt` contains a series of coordinates serialized as pairs of axial and radial coordinates of the cod-end meridian (one line of twines). Multiple runs are appended to the same result files and each column corresponds to one simulation from `results.txt`. Since the structure us axis-symmetric this is sufficient to recover the final deformed sape of the cod-end. In `shapes.txt`. As a sanity check notice that the 2nd coordinate (radial coordinate of the first node) in both columns equals to the entrance radius. The second to last coordinate (axial coordinate of the last node) equals to the result cod-end length.

```
0.00000E+000   0.00000E+000   
4.00000E-001   4.00000E-001   
1.30940E-002   1.31977E-002   
3.98180E-001   4.00862E-001   
...
1.16777E+001   1.04468E+001   
0.00000E+000   0.00000E+000
```

## Dependencies
- [Accord.Math]() - dense matrix computations (will be removed soon)
- [CSparse]() - sparse matrix computations
- [YamlDotNet]() - yaml job configuration parsing
- [CommandLineParser]() - CLI argments parsing