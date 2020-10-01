# Import OVITO modules.
import ovito
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
from ovito.vis import *

# Import Qt modules
from PySide2.QtGui import QPainter
from PySide2.QtGui import QFont
from PySide2.QtGui import QColor

from timeit import default_timer as timer

# Import standard Python and NumPy modules.
import sys
import numpy as np
import glob
import os 
from pathlib import Path

print("Opend file %s..." % Path(__file__).absolute())

print("ovito: %i.%i.%i" % ovito.version)
print("NumPy: {}".format(np.__version__))

#=== Load input data and create an ObjectNode with a data pipeline.

filepath = 'C:/Users/USER/Desktop/BionicMaterial/script/dump/'
atomdump = 'dump.Si.*'
bonddump = 'bond.dump.*'

firstn = 99

print("Reading atoms dump file...")

pipeline = import_file(filepath+atomdump, columns =
		 ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z", "Potential Energy", "StressTensor.XX", "StressTensor.YY", "StressTensor.XY"], multiple_frames = False, sort_particles=True)

print("StressTensor.XX", "StressTensor.YY", "StressTensor.XY")
#import stress= "Stress Tensor.XX", "Stress Tensor.YY", "Volumetric Stress"

#=== Display settings

atomtypes = pipeline.source.data.particles.particle_types

atomtypes.type_by_id(1).name = 'A' # stiff
atomtypes.type_by_id(1).radius = 75
atomtypes.type_by_id(1).color = np.array([255, 0, 0])/255.0

#atom沒有定義tpye 2
#atomtypes.type_by_id(1).name = 'B' # soft
#atomtypes.type_by_id(1).radius = 75
#atomtypes.type_by_id(1).color = np.array([0, 0, 255])/255.0

#=== Import bond information
bonds_dict = {}



print("Reading bonds dump file...")

start = timer()
atoms_files = sorted(glob.glob(filepath+atomdump), key = os.path.getmtime)
#print(atoms_files)
bonds_files = sorted(glob.glob(filepath+bonddump), key = os.path.getmtime)
#print(bonds_files)
for i in range(len(atoms_files)):
	t = pipeline.compute(frame=i).attributes['Timestep']
	print('\tReading frame: {:5d}/{:5d}'.format(i,len(atoms_files)),end="\r")
	print(i, bonds_files[i])
	with open(bonds_files[i],'r') as bonds_file:
		lines = bonds_file.readlines()
		t = 0
		n_bonds = 0
		for n, line in enumerate(lines):
			if("TIMESTEP" in line):
				t = int(lines[n+1])
			elif("NUMBER OF ENTRIES" in line):
				nbonds = int(lines[n+1])
			elif("ITEM: ENTRIES" in line):
				bonds = [ list(map(float, line.strip("\n").split()[1:5] )) 
             for line in lines[(n+1):(n+1+nbonds)]]
				bonds_dict[t] = np.reshape(bonds, (len(bonds),4))
		del lines

end = timer()
print("Elapsed Time: ", end-start)

#CUSTOM Modifier that creates bonds object for every frame
def create_bonds_modifier(frame, data):
	#find closest matching timestep:
	t = int(data.attributes["Timestep"])
	if( t not in bonds_dict.keys()):
		t = min(bonds_dict.keys(), key=lambda x:abs(x-t))
	data.particles_.bonds = Bonds()
	bonds = bonds_dict[t]
	data.particles_.bonds_.create_property('Topology', data=bonds[:,1:2]-np.ones((bonds.shape[0],2)))
	data.particles_.bonds_.create_property('Bond Type', data=bonds[:,2])
	data.particles.bonds.vis.width = 60
	data.particles.bonds.vis.shading = BondsVis.Shading.Flat

pipeline.modifiers.append(create_bonds_modifier)


pipeline.modifiers.append(ComputePropertyModifier(
    output_property = 'Volumetric Stress',
    expressions = ['(StressTensor.XX + StressTensor.YY + StressTensor.XY)/3.']
))


pipeline.modifiers.append(AtomicStrainModifier(
	cutoff = 200,
	output_strain_tensors = True
))

freeze = FreezePropertyModifier(source_property = 'Position',
								  destination_property = 'Position',
								  freeze_at = 0)
#freeze property 是產生裂縫的

pipeline.modifiers.append(freeze)

pipeline.add_to_scene()

vis_particle = pipeline.source.data.particles.vis
vis_particle.enabled = False

#vis_bond = pipeline.source.data.particles.bonds.vis
#vis_bond.enabled = True

vis_cell = pipeline.source.data.cell.vis
vis_cell.enabled = False


vp = Viewport(type = Viewport.Type.Top,
			  camera_pos = (12800,12800,0.5),
			  fov = 12800)

vp.render_anim(filename=filepath+'data/crack.png',
			   size=(128, 128),
			   background=(0, 1, 0),
			   range=(0,firstn))

for frame_index in range(0,pipeline.source.num_frames,2):
    index = int(frame_index/2)
    vp.render_image(filename=filepath+'data/crack{:04d}.png'.format(index),
                    size=(128,128),frame=frame_index,background=(0, 1, 0))


pipeline.modifiers.append(ColorCodingModifier(
    property = 'Volumetric Stress',
    gradient = ColorCodingModifier.BlueWhiteRed(),
	start_value = -3e3,
	end_value = 3e3
))


vis_particle.enabled = True

pipeline.compute()

vp.render_anim(filename=filepath+'data/stress.png',
			   size=(128, 128),
			   background=(0, 0, 0),
			   range=(0,firstn))

for frame_index in range(0,pipeline.source.num_frames,2):
    index = int(frame_index/2)
    vp.render_image(filename=filepath+'data/stress{:04d}.png'.format(index),
                    size=(128,128),frame=frame_index,background=(0, 0, 0))

pipeline.modifiers.append(ColorCodingModifier(
    property = 'Volumetric Strain',
    gradient = ColorCodingModifier.BlueWhiteRed(),
	start_value = -1,
	end_value = 1
))


vis_particle.enabled = True

pipeline.compute()

vp = Viewport(type = Viewport.Type.Top,
			  camera_pos = (12800,12800,0.5),
			  fov = 12800)

vp.render_anim(filename=filepath+'data/strain.png',
			   size=(128, 128),
			   background=(0, 0, 0),
			   range=(0,firstn))

for frame_index in range(0,pipeline.source.num_frames,2):
    index = int(frame_index/2)
    vp.render_image(filename=filepath+'data/strain{:04d}.png'.format(index),
                    size=(128,128),frame=frame_index,background=(0, 0, 0))
