import math
import numpy as np
import os
from simulationSettingsMicro import Names



'''
    ABAQUS part
'''



import abaqus
import abaqusConstants as const
import odbAccess


job_names = []
for mdl_name in Names.MODEL_NAMES:
    job_name = Names.JOB_NAME + '_' + mdl_name
    job_names.append(job_name)

Names.JOB_NAMES = job_names

def calc_avgStress(frame,coordSys,frameIVOL,totVolume,strainMeasure):
    '''
        It takes the name of the step and will return the average stress and strain
        args:
                step_name:str
        returns:
                Avg_Stress, Avg_Strain
    '''
    stress_transformed = frame.fieldOutputs['S'].getTransformedField(datumCsys=coordSys) 
    strain_transformed = frame.fieldOutputs[strainMeasure].getTransformedField(datumCsys=coordSys) 
    stress_transformed = stress_transformed.getSubset(position=const.INTEGRATION_POINT)
    strain_transformed = strain_transformed.getSubset(position=const.INTEGRATION_POINT)
    # Stress Sum
    Tot_Stress = np.sum([stress_transformed.values[II].data * frameIVOL.values[II].data for II in range(len(stress_transformed.values))],0)
    # Strain Sum
    Tot_Strain = np.sum([strain_transformed.values[II].data * frameIVOL.values[II].data for II in range(len(strain_transformed.values))],0)
    # Calculate Average
    Avg_Stress = Tot_Stress/totVolume
    Avg_Strain = Tot_Strain/totVolume
    print('Average stresses Global CSYS: 11-22-33-12-13-23')
    # print('Avg_Strain: ',Avg_Strain)
    return Avg_Stress, Avg_Strain

def writeDataAsCSV():
    file_name = Names.STRESS_STRAIN_DATA_FILE_NAME + '.csv'
    with open(file_name,'w') as f:
        for i, job_name in enumerate(Names.JOB_NAMES):
            job_name = job_name + '.odb'
            odb = odbAccess.openOdb(path=job_name, readOnly=True, readInternalSets=True)
            stiffnessHomogStep = odb.steps[odb.steps.keys()[0]]
            all_frames = [ stiffnessHomogStep.frames[i] for i in range(1, len(stiffnessHomogStep.frames)) ]
            frameIVOL = stiffnessHomogStep.frames[0].fieldOutputs['IVOL'].getSubset(
                                                position=const.INTEGRATION_POINT)
            totVolume = sum([frameIVOL.values[i].data for i in range(len(frameIVOL.values))])
            print('IPVolume = %1.4f'%totVolume)
            
            actual_volume = 8*Names.a1 * Names.a2 * Names.a3
            print('Actual vol7ume = %1.4f'%actual_volume)
            print('abs diff vol: %1.4e' %(actual_volume- totVolume))
            coordSysName = 'global'
            coordSys = odb.rootAssembly.DatumCsysByThreePoints(name=coordSysName, 
                                                            coordSysType=const.CARTESIAN, 
                                                            origin=(0.0,0.0,0.0), 
                                                            point1=(1.0,0.0,0.0), 
                                                            point2=(0.0,1.0,0.0))
            preferredStrainMeasures = ('NE', 'E', 'LE')
            # preferredStrainMeasures = ('LE','E')
            strainMeasure = None
            for measure in preferredStrainMeasures:
                if measure in stiffnessHomogStep.frames[1].fieldOutputs:
                    strainMeasure = measure
                    break
            print('prefered strain measure for job=%s is %s'%(job_name,strainMeasure))
            for step_name in odb.steps.keys():
                stiffnessHomogStep = odb.steps[step_name]
                all_frames = [ stiffnessHomogStep.frames[i] for i in range(1, len(stiffnessHomogStep.frames)) ]
                frameIVOL = stiffnessHomogStep.frames[0].fieldOutputs['IVOL'].getSubset(
                                                position=const.INTEGRATION_POINT)
                for frame_num, frame in enumerate(all_frames):
                    try:
                        avg_stress, avg_strain = calc_avgStress(frame=frame,
                                                                coordSys=coordSys,
                                                                frameIVOL=frameIVOL,
                                                                totVolume=totVolume,
                                                                strainMeasure = strainMeasure)
                        for strain_element in avg_strain:
                            f.write(str(strain_element) + ', ' )
                        for stress_element in avg_stress:
                            f.write(str(stress_element) + ', ' )
                        f.write('\n')
                        print('frame=%i is completed!'%(frame_num+1))
                    except:
                        print('There was error in frame %i and in step %s'%(frame_num,step_name))
                print('step=%s is completed!'%(step_name))
            

writeDataAsCSV()