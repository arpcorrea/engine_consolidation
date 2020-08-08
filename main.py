from consolidate import *



cwd = os.getcwd()
# All the data from deck.yaml is now in the following deck variable

deck = Deck( cwd + "/AL_problem.yaml" )

geometry = Geometry(deck)

meshes = Plate( deck,geometry )

model_HT= HeatTransfer(deck,meshes)

# model_IC = IntimateContact(meshes,deck)

plots=PlotsTwoPlates(deck,meshes,meshes.T)

solves = SolvesTwoPlates( deck,model_HT,meshes,plots)
