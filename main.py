from consolidate import *



cwd = os.getcwd()
# All the data from deck.yaml is now in the following deck variable

deck = Deck( cwd + "/AL_problem.yaml" )

geometry = Geometry(deck)

meshes = MeshOnePlate( deck,geometry )

BC = BoundaryCondition(deck, geometry, meshes)

model_HT= HeatTransfer(geometry,deck,meshes, BC)

plots=PlotsTwoPlates(deck,meshes,BC)

solves = SolvesOnePlates( deck, meshes, BC,model_HT,plots)
