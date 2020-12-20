from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

pdb = PDBFile('initial.pdb')
forcefield = ForceField('ff.xml')
system = forcefield.createSystem(pdb.topology)

#Agregar fuerzas
force=CustomNonbondedForce('energy;'
                           'energy=4*epsilon*((sigma/r)^12-(sigma/r)^6);'
                           'sigma=0.5*(sigma1+sigma2);' 
                           'epsilon=sqrt(epsilon1*epsilon2);')

force.addPerParticleParameter('sigma')
force.addPerParticleParameter('epsilon')
epsilon=-1
sigma=1
for i in range(50):
    force.addParticle([epsilon,sigma])
system.addForce(force)

d1=4.101428012677189
d2=5.685486977714755
d3=8.912940586876381
k=100

force= HarmonicBondForce()
for i in range(49):
  force.addBond(i,i+1,d1,k)
for i in range(48):
  force.addBond(i,i+2,d2,k)
for i in range(47):
  force.addBond(i,i+3,d3,k)
system.addForce(force)

integrator = LangevinIntegrator(300*kelvin, 1E-4/picosecond, 0.1*picoseconds)
platform = Platform.getPlatformByName('OpenCL') #Cambiar a CPU o CUDA dependiendo de donde esten corriendo
simulation = Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('LJ_output.pdb', 100000))
simulation.reporters.append(StateDataReporter(stdout, 100000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(100000000)
