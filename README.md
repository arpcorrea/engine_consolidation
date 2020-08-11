# The problem


A plate with dimensions X and Y are defined by the user in ["Geometry"]["Length X"]["Length Y"]

The plate, made of orthotropic, material is heated on the left boundary with a constant heat Q, defined by the user in ["Boundary Condition"]["Input Power Density"]. We assumed that this material is aluminium"

The plate's initial temperature is defined by the user in ["Boundary Condition"]["Initial Temperature"]. The material properties are defined in ["Materials"]["Aluminium"].

The plate is inside a room, with a constant temperature is set by the user in ["Boundary Condition"]["Room Temperature"]. The convection coefficient is defined in ["Boundary Condition"]["Convective Coefficient"]

The user also define the number of elements in X and Y in ["simulation"]["Number of Elements X"] and ["simulation"]["Number of Elements Y"] respectivelly.

The number of steps and step time are defined in ["simulation"]["Step Time"] and ["simulation"]["Number Time Steps"] respectivelly. Notice that large step time may not converge and an error message will popup.

Finally, the code generates Temperature plots at every interval, defined in ["Plot"]["plot interval"].



```yaml
Problem Type:
    Type: Aluminium
    
Geometry:
    Length X: 0.2
    Length Y: 0.2
 
Materials:
  Aluminium:
    Thermal Conductivity X: 237
    Thermal Conductivity Y: 237
    Density: 2700
    Specific Heat Capacity: 1000
    
Boundary Condition:
    Initial Temperature: 293
    Input Power Density: 2613200
    Room Temperature: 293
    Convective Coefficient: 150
    Ideal Temperature: 1200

Simulation:
  Time Step: 0.01
  Number Time Steps: 20001
  Number of Elements X: 100 
  Number of Elements Y: 100

Plot:
  Temp Output Folder: "./output/Temperature/"
  figure temperature name: Temperature
  Color Interpolation: 50
  Color Map: "inferno"
  plot interval: 500
  
Animation:
    name: temperature
```


# Getting Started


### Geomtry
A class named "Geometry" is difined only to store the X and Y dimensions of the problem.
### Mesh
A class named "MeshOnePlate" is defined to create the mesh.
### BC
a class named "Boundary Condition" assigns Initial Temperature, Thermal Diffusivity in X and Thermal Diffusivity in Y, Cp, Density in each element. 
### Model
A class named "Heat Transfer" is the thermal engine. It solves convection on the edges and the heat transfer problem by central differences within the plate.
### Plot and Animation
A class named "PlotsTwoPlates" is defined to generate figures of the thermal history in spaced intervals. In the end, a .GIF is created from those figures.
### Solver
A class named "SolvesTwoPlates" is defined to create a forward step analysis, solving the heat transfer model and plotting the figures when convenient.