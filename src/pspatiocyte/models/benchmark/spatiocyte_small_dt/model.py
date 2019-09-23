import time
import math
import numpy as np

nx = 94
ny = 94
nz = 94
D = 10e-12 # m^2/s
volume = 909e-18 # m^3
rv = math.pow(volume/((nz*2.0)*(ny*math.sqrt(3.0))*(nx*math.sqrt(8.0/3.0))),
        1.0/3.0)
box_x = rv*(nx-3)*math.sqrt(8.0/3.0)
box_y = rv*(ny-3)*math.sqrt(3.0)
box_z = rv*(nz-3)*2.0

sim = theSimulator.createStepper('SpatiocyteStepper', 'SS')
sim.VoxelRadius = rv;

theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:GEOMETRY').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = box_x;
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = box_y;
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = box_z;
#periodic boundaries
theSimulator.createEntity('Variable', 'Variable:/:XYPLANE').Value = 1;
theSimulator.createEntity('Variable', 'Variable:/:XZPLANE').Value = 1;
theSimulator.createEntity('Variable', 'Variable:/:YZPLANE').Value = 1;

theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:E').Value = 9090
theSimulator.createEntity('Variable', 'Variable:/:S').Value = 90910
theSimulator.createEntity('Variable', 'Variable:/:ES').Value = 0
theSimulator.createEntity('Variable', 'Variable:/:P').Value = 0

#log = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:log')
#log.VariableReferenceList = [['_', 'Variable:.:A']]

pop = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:pop')
pop.VariableReferenceList = [['_', 'Variable:.:E']]
pop.VariableReferenceList = [['_', 'Variable:.:S']]

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffE')
dif.VariableReferenceList = [['_', 'Variable:.:E']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffS')
dif.VariableReferenceList = [['_', 'Variable:.:S']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffES')
dif.VariableReferenceList = [['_', 'Variable:.:ES']]
dif.D = D

dif = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffP')
dif.VariableReferenceList = [['_', 'Variable:.:P']]
dif.D = D

keff = 0.01e-18
kd = 4.*np.pi*2*rv*2*D
ka = keff*kd/(kd-keff)
binder = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:fwd')
binder.VariableReferenceList = [['_', 'Variable:.:E','-1']]
binder.VariableReferenceList = [['_', 'Variable:.:S','-1']]
binder.VariableReferenceList = [['_', 'Variable:.:ES','1']]
binder.k = ka

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:back')
react.VariableReferenceList = [['_', 'Variable:.:ES', '-1']]
react.VariableReferenceList = [['_', 'Variable:.:E', '1']]
react.VariableReferenceList = [['_', 'Variable:.:S', '1']]
react.k = 1

react = theSimulator.createEntity('SpatiocyteNextReactionProcess', 'Process:/:prod')
react.VariableReferenceList = [['_', 'Variable:.:ES', '-1']]
react.VariableReferenceList = [['_', 'Variable:.:E', '1']]
react.VariableReferenceList = [['_', 'Variable:.:P', '1']]
react.k = 1

log = theSimulator.createEntity('IteratingLogProcess', 'Process:/:log')
log.VariableReferenceList = [['_', 'Variable:.:E']]
log.VariableReferenceList = [['_', 'Variable:.:S']]
log.VariableReferenceList = [['_', 'Variable:.:ES']]
log.VariableReferenceList = [['_', 'Variable:.:P']]
log.LogInterval = 0.01
log.LogEnd = 10
log.FileName = "output.txt"

start = time.time()
run(10.0)
end = time.time()
print "total runtime:", end-start
