import math
import os 
import sys
import numpy as np

os.chdir(r'%s'%os.getcwd())
###################################
#######  Calculating the geometrical parameters
# num_fibers = 1

# radius of fiber
unit_conversion_factor = 1.0
r_f = 300.0e-3 * unit_conversion_factor
# volume fraction
volume_f = 0.4
# spacing between chords 
spacing_chords = 0.738117 * unit_conversion_factor

# volume of fiber = pi * (r_f)**2 * 2*a1
# volume of fibers = num_fibers * volume of fiber = 4 * a1 * pi * (d_f/2)**2
# volume of RVE = 2a1 * 2a2 * 2a3 = 8 * a1 * a2 * a3


a3 = spacing_chords/2.0

a2 = math.pi/(4.0*a3 * volume_f) * (r_f)**2

# a1 can be chosen as arbitarely 
a1 = a2/4.0
''' Mesh size parameter change it as per your analysis'''
mesh_size = spacing_chords/40.0

print('a1 = %f, (length of fiber)                   2 a1=%f'%(round(a1,3),round(2*a1,3)))
print('a2 = %f, (vertical spacing between fibers)   2 a2=%f'%(round(a2,3),round(2*a2,3)))
print('a3 = %f, (horizontal spacing between fibers) 2 a3=%f'%(round(a3,3),round(2*a3,3)))


# spcifying the other point for the fiber

'''
    ABAQUS part
'''
from simulationSettingsMicro import Names
from mechanicalLoading import pipline_3D_periodicity
from materialLibrary import create_materials
from initialSingleChordGeometry import pipeLineInitialGeometry
Names.a1 = a1
Names.a2 = a2
Names.a3 = a3
Names.r_f = r_f


model_name    = Names.MODEL_NAME
part_name     = Names.PART_NAME
instance_name = Names.INSTANCE_NAME
load_case_name     = Names.STEP_NAME

import abaqus
import abaqusConstants as const
import part
import assembly
import interaction
import mesh


# from abaqusConstants import *

# Creating a cube model
# if 'Model-1' exists then rename it to eps11
if 'Model-1' in abaqus.mdb.models.keys():
    abaqus.mdb.models.changeKey(fromName=model_name, toName='eps11')
    Names.MODEL_NAME = 'eps11'
    model_name = Names.MODEL_NAME


eps11_model = abaqus.mdb.models[model_name]
abaqus.session.journalOptions.setValues(replayGeometry=const.COORDINATE,recoverGeometry=const.COORDINATE)
# creating initial geometry
pipeLineInitialGeometry(eps11_model,a1,a2,a3,r_f)


# create materials
create_materials(createdModel=eps11_model)
material_names = {'fiber':'GEW-661344-MARLOW', 'matrix':'BONN-60SHA'}

# creating sections
eps11_model.HomogeneousSolidSection(material=material_names['fiber'], name='Fiber', 
    thickness=None)
eps11_model.HomogeneousSolidSection(material=material_names['matrix'], name='Matrix', 
    thickness=None)
# Creating a set for the fiber and matrix
p = eps11_model.parts[part_name]
p.Set(cells=
    p.cells.findAt(
        ((0.0,0.0,a3),),  
        )
    , name=
    'Fibre')
p.Set(cells=
    p.cells.findAt(
        ((0.0,0.0,0.0),),  
        )
    , name=
    'Matrix')
# Asigning set 1 to the fiber section
p.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=const.MIDDLE_SURFACE, region=
    p.sets['Fibre'], sectionName='Fiber', 
    thicknessAssignment=const.FROM_SECTION)
p.SectionAssignment(offset=0.0, 
    offsetField='', offsetType=const.MIDDLE_SURFACE, region=
    p.sets['Matrix'], sectionName='Matrix', 
    thicknessAssignment=const.FROM_SECTION)

# creating assembly instant
eps11_model.rootAssembly.DatumCsysByDefault(const.CARTESIAN)
eps11_model.rootAssembly.Instance(dependent=const.OFF, name=instance_name, part=p)

# step
def createFiniteSteps(model,num_steps):
    # model.StaticLinearPerturbationStep(name= 'Step-1', previous='Initial')
    model.StaticStep(name='Step-1', previous='Initial', nlgeom=const.ON,
                     maxNumInc=10000, 
                     initialInc=0.001,
                     minInc = 1e-20)
    # channging default time incrmination from 5 to 25
    model.steps['Step-1'].control.setValues(timeIncrementation=(4.0, 
        8.0, 9.0, 16.0, 10.0, 4.0, 12.0, 25.0, 6.0, 3.0, 50.0),resetDefaultValues=const.OFF,allowPropagation=const.OFF)
    for i in range(num_steps-1):
        # model.StaticLinearPerturbationStep(name='Step-'+str(i+2), previous='Step-'+str(i+1))
        current_step_name = 'Step-'+str(i+2)
        model.StaticStep(name= current_step_name, previous='Step-'+str(i+1),
                         maxNumInc=10000, 
                         initialInc=0.001,
                         minInc = 1e-20)
        model.steps[current_step_name].control.setValues(timeIncrementation=(4.0, 
            8.0, 9.0, 16.0, 10.0, 4.0, 12.0, 25.0, 6.0, 3.0, 50.0),resetDefaultValues=const.OFF,allowPropagation=const.OFF)



from geometryCreationSingleChord import partitionGeometry
partitionGeometry(Model=eps11_model,r_f=r_f,a1=a1,a2=a2,a3=a3,Names=Names)


# Meshing
def meshByQuad(createdModel=eps11_model,instance_name=instance_name,mesh_size=mesh_size):
    elemType1 = mesh.ElemType(elemCode=const.C3D8RH, elemLibrary=const.STANDARD, 
        kinematicSplit=const.AVERAGE_STRAIN, hourglassControl =const.DEFAULT)
    elemType2 = mesh.ElemType(elemCode=const.C3D6, elemLibrary=const.STANDARD)
    elemType3 = mesh.ElemType(elemCode=const.C3D4, elemLibrary=const.STANDARD)
    a = createdModel.rootAssembly
    c1 = a.instances[instance_name].cells
    a.Set(cells= c1,name='all-cells') 
    a.setElementType(regions=a.sets['all-cells'], elemTypes=(elemType1, elemType2, 
        elemType3))
    partInstances =(a.instances[instance_name], )
    a.seedPartInstance(regions=partInstances, size=mesh_size, deviationFactor=0.1, 
        minSizeFactor=0.1)
    a.generateMesh(regions=partInstances)

# createMeshBySweepOneChord(createdModel=createdModel,instance_name=instance_name,mesh_size=mesh_size)
meshByQuad(createdModel=eps11_model,instance_name=instance_name,mesh_size=mesh_size)


created_models = {}
created_models[Names.MODEL_NAMES[0]]=eps11_model
# copy eps11 model to create other models
for i in range(1,6):
    created_models[Names.MODEL_NAMES[i]]=abaqus.mdb.Model(name=Names.MODEL_NAMES[i], objectToCopy=abaqus.mdb.models['eps11'])

Names.LAST_STEP_NAMES = [1,1,1,1,1,1] 
for i in range(6):
    createFiniteSteps(model=created_models[Names.MODEL_NAMES[i]],num_steps=Names.NUM_STEPS_EACH_MODEL[i])
    # output requests
    created_models[Names.MODEL_NAMES[i]].fieldOutputRequests['F-Output-1'].setValues(variables=(
        'NE','LE','E', 'S', 'U', 'IVOL'))
    Names.LAST_STEP_NAMES = created_models[Names.MODEL_NAMES[i]].steps.keys()[-1]
    pipline_3D_periodicity(Model=created_models[Names.MODEL_NAMES[i]],Names=Names,model_idx=i)

Names.LAST_STEP_NAME = eps11_model.steps.keys()[-1]




# create job
# create job for all models
job_names = []
for mdl_name in Names.MODEL_NAMES:
    job_name = Names.JOB_NAME + '_' + mdl_name
    job_names.append(job_name)
    abaqus.mdb.Job(atTime=None, contactPrint=const.OFF, description='', echoPrint=const.OFF, 
            explicitPrecision=const.SINGLE, getMemoryFromAnalysis=True, historyPrint=const.OFF, 
            memory=90, memoryUnits=const.PERCENTAGE, model=mdl_name, modelPrint=const.OFF, 
            multiprocessingMode=const.DEFAULT, name=job_name, nodalOutputPrecision=const.SINGLE, 
            numCpus=1, numGPUs=0, queue=None, resultsFormat=const.ODB, scratch='', type=
            const.ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

Names.JOB_NAMES = job_names
# submit job
for job_name in job_names:
    abaqus.mdb.jobs[job_name].submit(consistencyChecking=const.OFF)

# waiting for it to complete
for job_name in job_names:
    abaqus.mdb.jobs[job_name].waitForCompletion()
# create PINNED BC at the center of the RVE Or do not create it
execPyFile('writeStressStrainData.py')
