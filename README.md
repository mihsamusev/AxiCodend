# How to build
Install the latest .NET Core SDK
Install Git
Clone this repo
`build.sh` in the root of the cloned repo

# How i migrated
```sh
mkdir axi_codend
dotnet new console
# copy all files and overwrite Program.cs

# add dependencies
dotnet add package Accord.Math --version 3.8.2-alpha # is obsolete, will be refactored later
dotnet add package CSparse --version 3.6.0
dotnet add package YamlDotNet --version 11.2.1

dotnet build
```

# Main refractoring ideas
- MeshOrientation is an enum
- inputs are JSON
- Standart solver settings are JSON
- outputs are JSON

```C#
var cod = new Codend()
    .type()
    .material()
    .build()

var sim = new Simulation()
    .with_codend(cod)
    .with_towing()
    .with_catch()
    .build()
```

```sh
# simulate towing x catch matrix
axi_codend job.json
axi_codend job.json -s solver_settings.json
```

```C#
// parse job.json args

// build codend (builder is independent on parser)
var codend = CodendBuilder().FromJson(filename).build()
// or
var codend = CodendBuilder()
    .withMeshesAlong(l)
    .withMeshesAcross(nr)
    .withRadius(r)
    .WithMaterial(material).build()

var sim = SimulationBuilder()
    .withTowing(towing)
    .withCatch(catch)
    .withSettings(solverSettings)

sim.simulate()

sim.save_json(path)
```
# Redundant fucntions

UpdatePosition()
ApplyTowing()
UpdateTotalForces()
UpdateResidual()
UpdateTotalJacobian()

# AxiCodend
Simulation of trawl cod-end shapes during the towing using an object-oriented C# version of axis-symmetric numerical model developed by D.Priour (see atached paper and thesis). Input happens throught the "input.txt" file, and should follow the format:

*Block of paths (essential input) contains the paths where the results are stored, the input is by default input.txt and has to be in the same directory as .exe file*
```
Paths
OutputShapes  C:\Users\MyFavouriteFolder\filewithshapes.txt
OutputResults C:\Users\MyFavouriteFolder\filewithresults.txt
```
*Block of materials (essential input)*
```
Material       
MeshSide            0.120000            
KnotSize            0.012000            
TwineEA             1000.000000         
KnotEA              1000.000000         
MeshOrientation     90                  
```
*Block of cod-end structure (essential input)*
```
Codend         
MeshesAlong         100                 
MeshesAround        100                 
EntranceRadius      0.400000            
```
*Block of catches (essential input), represents how many meshes along are blocked by catch. Make sure that Count matches the amount of elements in the following catches column.*
```
Catch          
Count               2                  
50
75
```
Block of Towing (essential input)
```
Towing         
TowingSpeed         1.500000            
```
*Block of Solver settings (optional), if nothing is input the solver uses following default settings*
```
Solver         
IterMax         5000
ResidualTol     1e-3                          
DisplacementTol 1e-4                      
ResidualMax     1e20              
StiffnessTol    1            
DiagStiffness   1          
ReduceStiffnessBy   0.1
IncreaseStiffnessBy 2
MinCatchBlock       5
MinTowingSpeed      0.1
UsePreviousAsPrecalc 1
```
