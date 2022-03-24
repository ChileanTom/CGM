from __future__ import print_function
###################################################################
readParamsFromTable(rParticle=0.09, rRelFuzz=0,rCoff=2,bot_limit=1,width=1,tot_limit=1,toy_limit=1,ejex_limit=5.5,num_spheres=500,thick = 0.01,key='_define_a_name_',stabilityThreshold=0.5)
from yade.params.table import *
from numpy import arange
from yade import pack
import pylab
from yade import pack
import gts, os.path, locale
# -*- encoding=utf-8 -*-
""" CAUTION:
Running this script can take very long!
"""
locale.setlocale(
        locale.LC_ALL, 'en_US.UTF-8'
)  #gts is locale-dependend.  If, for example, german locale is used, gts.read()-function does not import floats normally
'''
if you get "Error: unsupported locale setting"
-> type as root: "dpkg-reconfigure locales"
-> choose "en_US.UTF-8" (press space to choose)
'''

############################################
### DEFINING VARIABLES AND MATERIALS ###
############################################
## corners of the initial packing
mn,mx=Vector3(0,0,0),Vector3(5.5,5.5,1)
compFricDegree = 30
finalFricDegree = 30
rate=0.2
damp=0.1
young=50e6 # contact stiffness
## create material #0, which will be used as default
#O.materials.append(FrictMat(young=5e6,poisson=0.5,frictionAngle=radians(compFricDegree),density=2600,label='spheres'))
O.materials.append(FrictMat(young=5e6,poisson=0.5,frictionAngle=0,density=0,label='walls'))
## create walls around the packing
walls=aabbWalls([mn,mx],thickness=thick,material='walls')
wallIds=O.bodies.append(walls)
#Define Materials
Rockfill=O.materials.append(FrictMat(young=100e6,poisson=0.3,density=2650,frictionAngle=radians(30),label='spheres'))
# define the section shape as polygon in 2d; repeat first point at the end to close the polygo
surf = gts.read(open('talud2.coarse.gts'))
# fill this solid with triaxial packing; it will compute minimum-volume oriented bounding box
# to minimize the number of throw-away spheres.
# It does away with about 3k spheres for rParticle 3e-2
sp1 = SpherePack()
sp1 = pack.randomDensePack(pack.inGtsSurface(surf), radius=rParticle,material=Rockfill, rRelFuzz=rRelFuzz,spheresInCell=1000, memoizeDb='/tmp/gts-triax.sqlite', returnSpherePack=True)
rockfill = sp1.toSimulation()

######################################################################################################3
#################################################################################
#AXIS Y (vertical)
#AXIS X (up-down stream)
#AXIS Z (width)
bot = [O.bodies[s] for s in rockfill if O.bodies[s].state.pos[1]<rParticle*rCoff*2]
tot = [O.bodies[s] for s in rockfill if O.bodies[s].state.pos[2]<=rParticle*rCoff*2]
toy = [O.bodies[s] for s in rockfill if O.bodies[s].state.pos[2]>=width-rParticle*rCoff*2]
ejex = [O.bodies[s] for s in rockfill if O.bodies[s].state.pos[0]<=rParticle*rCoff*2]
ejexx = [O.bodies[s] for s in rockfill if O.bodies[s].state.pos[0]>=ejex_limit-rParticle*rCoff*2]

for s in bot:
        if s.state.pos[1]<=bot_limit:
                bot_limit = s.state.pos[1]#Define the minimal position Y (vertical) from the dense particles
                bot_id = s.id
for b in bot: 
        b.state.blockedDOFs = 'xyzXYZ'
        b.state.vel = (0,0,0)

for s in tot:
        if s.state.pos[2]<=tot_limit:
                tot_limit = s.state.pos[2]
                tot_id = s.id
for b in tot: #reemplazar wall por top o tot layers
        b.state.blockedDOFs = 'xyzXYZ'
        b.state.vel = (0,0,0)

for s in toy:
        if s.state.pos[2]<=toy_limit:
                toy_limit = s.state.pos[2]
                toy_id = s.id
for b in toy: #reemplazar wall por top o toy layers
        b.state.blockedDOFs = 'xyzXYZ'
        b.state.vel = (0,0,0)

for s in ejex:
        if s.state.pos[0]<=ejex_limit:
                ejex_limit = s.state.pos[0]
                ejex_id = s.id
for b in ejex: #reemplazar wall por top o ejex layers
        b.state.blockedDOFs = 'xyz'
        b.state.vel = (0,0,0)

for s in ejexx:
        if s.state.pos[0]<=ejex_limit:
                ejexx_limit = s.state.pos[0]
                ejexx_id = s.id
for b in ejexx: #reemplazar wall por top o ejexx layers
        b.state.blockedDOFs = 'xyz'
        b.state.vel = (0,0,0)

O.dt=.5*PWaveTimeStep() # initial timestep, to not explode right away
O.usesTimeStepper=True

############################
### DEFINING ENGINES ###
############################
triax=ThreeDTriaxialEngine(
	maxMultiplier=1.+2e4/young, # spheres growing factor (fast growth)
	finalMaxMultiplier=1.+2e3/young, # spheres growing factor (slow growth)
	thickness = thick,
	stressControl_1 = True,
	stressControl_2 = True,
	stressControl_3 = True,
	## Independant stress values for anisotropic loadings
	goal1=-10000,
	goal2=-10000,
	goal3=-10000,
	internalCompaction=True,
	Key=key,
)

newton=NewtonIntegrator(damping=damp)

O.engines=[
 ForceResetter(),
 InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
 InteractionLoop(
  [Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
  [Ip2_FrictMat_FrictMat_FrictPhys()],
  [Law2_ScGeom_FrictPhys_CundallStrack()]
 ),
 ## We will use the global stiffness of each body to determine an optimal timestep (see https://yade-dem.org/w/images/1/1b/Chareyre&Villard2005_licensed.pdf)
 GlobalStiffnessTimeStepper(active=1,timeStepUpdateInterval=100,timestepSafetyCoefficient=0.8),
 triax,
 TriaxialStateRecorder(iterPeriod=100,file='WallStresses'+key),
 newton
]


#Display spheres with 2 colors for seeing rotations better
Gl1_Sphere.stripes=0
yade.qt.Controller(), yade.qt.View()

#######################################
### APPLYING CONFINING PRESSURE ###
#######################################

while 1:
  O.run(1000, True)
  #the global unbalanced force on dynamic bodies, thus excluding boundaries, which are not at equilibrium
  unb=unbalancedForce()
  #average stress
  #note: triax.stress(k) returns a stress vector, so we need to keep only the normal component
  meanS=(triax.stress(triax.wall_right_id)[0]+triax.stress(triax.wall_top_id)[1]+triax.stress(triax.wall_front_id)[2])/3
  print('unbalanced force:',unb,' mean stress: ',meanS)
  if unb<stabilityThreshold and abs(meanS+10000)/10000<0.001:
    break

O.save('compressedState'+key+'.xml')
print("###      Isotropic state saved      ###")

##############################
### DEVIATORIC LOADING ###
##############################

#let us turn internal compaction off...
triax.internalCompaction=False

#
setContactFriction(radians(finalFricDegree))

#... and make stress control independant on each axis
triax.stressControl_1=triax.stressControl_2=triax.stressControl_3=True
# We have to turn all these flags true, else boundaries will be fixed
triax.wall_bottom_activated=True
triax.wall_top_activated=True#False
triax.wall_left_activated=True
triax.wall_right_activated=True
triax.wall_back_activated=True
triax.wall_front_activated=True


#If we want a triaxial loading at imposed strain rate, let's assign srain rate instead of stress
triax.stressControl_2=0 #we are tired of typing "True" and "False", we use implicit conversion from integer to boolean
triax.strainRate2=0.01
triax.strainRate1=triax.strainRate3=1000.0

O.run()

