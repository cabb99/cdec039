#Generacion de una configuracion aleatoria para el pdb
import pandas as pd
import numpy as np

n=100 #Numero de particulas
d=3 #Numero de dimensiones
box_A=1000 #Tamaño de la caja (A)

coords=np.random.random([n,d])*box_A
atoms=pd.DataFrame(coords,columns=['x','y','z'])
atoms['name']='A' #Nombre del atomo
atoms['resname']='AA' #Nombre del residuo
atoms['chainID']='A' #Nombre de la cadena
atoms['resSeq']=range(len(atoms)) #Numero del residuo
atoms['element']='O' #Elemento (tiene que ser un elemento verdadero)

pdb_file='one_particle.pdb'
# Write pdb file
with open(pdb_file, 'w+') as pdb:
    for i, atom in atoms.iterrows():
        pdb_line = f'ATOM  {i + 1:>5} {atom["name"]:^4} {atom.resname:<3} {atom.chainID}{atom.resSeq:>4}    {atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}' + ' ' * 22 + f'{atom.element:2}' + ' ' * 2
        assert len(pdb_line) == 80, 'An item in the atom table is longer than expected'
        pdb.write(pdb_line + '\n')

#Lennard-Jones gas
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

pdb = PDBFile('one_particle.pdb')
forcefield = ForceField('one_particle.xml')
system = forcefield.createSystem(pdb.topology)

#Tamaño de la caja (Sistema periodico)
box_nm=box_A/100
system.setDefaultPeriodicBoxVectors([box_nm,0,0],[0,box_nm,0],[0,0,box_nm])

#Agregar fuerzas
force=CustomNonbondedForce("""energy;
                              energy=4*epsilon*((sigma/r)^12-(sigma/r)^6); 
                              sigma=0.5*(sigma1+sigma2); 
                              epsilon=sqrt(epsilon1*epsilon2)""");

force.addPerParticleParameter('sigma')
force.addPerParticleParameter('epsilon')
epsilon=-10
sigma=10
for i in range(100):
    force.addParticle([epsilon,sigma])
system.addForce(force)

integrator = LangevinIntegrator(100*kelvin, 1/picosecond, 0.02*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('LJ_output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 10000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(1000000)
