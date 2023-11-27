
import abaqusConstants as const

def create_materials(createdModel):
    createdModel.Material(name='GEW-661344-MARLOW')
    createdModel.materials['GEW-661344-MARLOW'].Hyperelastic(
        materialType=const.ISOTROPIC, type=const.MARLOW, table=())
    createdModel.materials['GEW-661344-MARLOW'].hyperelastic.UniaxialTestData(
        smoothing=3, table=((0.0, 0.0), (18.066, 0.01), (32.699, 0.02), (
        46.648, 0.03), (61.936, 0.04), (79.935, 0.05), (101.442, 0.06), (
        126.756, 0.07), (155.751, 0.08), (187.953, 0.09), (222.615, 0.1), (
        258.796, 0.11), (295.432, 0.12), (331.412, 0.13), (365.659, 0.14), (
        397.198, 0.15), (425.237, 0.16), (449.241, 0.17), (469.009, 0.18), (
        484.746, 0.19), (500.695, 0.2033)))
    createdModel.materials['GEW-661344-MARLOW'].Density(table=((1.36e-09, 
        ), ))
    createdModel.Material(name='VLGQ6590BAY')
    createdModel.materials['VLGQ6590BAY'].Hyperelastic(
        materialType=const.ISOTROPIC, type=const.MARLOW, table=())
    createdModel.materials['VLGQ6590BAY'].hyperelastic.UniaxialTestData(
        smoothing=3, table=((0.0, 0.0), (0.0592828, 0.01), (0.115615, 0.02), (
        0.169171, 0.03), (0.220116, 0.04), (0.268607, 0.05), (0.314794, 0.06), 
        (0.358817, 0.07), (0.40081, 0.08), (0.440898, 0.09), (0.479201, 0.1), (
        0.51583, 0.11), (0.550892, 0.12), (0.584487, 0.13), (0.616708, 0.14), (
        0.647643, 0.15), (0.677375, 0.16), (0.705982, 0.17), (0.733537, 0.18), 
        (0.760108, 0.19), (0.785758, 0.2), (0.810548, 0.21), (0.834534, 0.22), 
        (0.857767, 0.23), (0.880296, 0.24), (0.902166, 0.25), (0.92342, 0.26), 
        (0.944096, 0.27), (0.964233, 0.28), (0.983862, 0.29), (1.00302, 0.3), (
        1.02173, 0.31), (1.04002, 0.32), (1.05792, 0.33), (1.07545, 0.34), (
        1.09264, 0.35), (1.10951, 0.36), (1.12607, 0.37), (1.14234, 0.38), (
        1.15834, 0.39), (1.1741, 0.4), (1.18961, 0.41), (1.20491, 0.42), (1.22, 
        0.43), (1.2349, 0.44), (1.24962, 0.45), (1.26417, 0.46), (1.27857, 
        0.47), (1.29283, 0.48), (1.30697, 0.49), (1.32099, 0.5), (1.3349, 
        0.51), (1.34873, 0.52), (1.36248, 0.53), (1.37616, 0.54), (1.3898, 
        0.55), (1.40339, 0.56), (1.41695, 0.57), (1.4305, 0.58), (1.44406, 
        0.59), (1.45762, 0.6), (1.47121, 0.61), (1.48483, 0.62), (1.49851, 
        0.63), (1.51226, 0.64), (1.52609, 0.65), (1.54001, 0.66), (1.55403, 
        0.67), (1.56818, 0.68), (1.58246, 0.69), (1.59688, 0.7), (1.61146, 
        0.71), (1.62622, 0.72), (1.64115, 0.73), (1.65628, 0.74), (1.67161, 
        0.75), (1.68716, 0.76), (1.70293, 0.77), (1.71892, 0.78), (1.73516, 
        0.79), (1.75163, 0.8), (1.76835, 0.81), (1.78532, 0.82), (1.80253, 
        0.83), (1.82, 0.84), (1.8377, 0.85), (1.85565, 0.86), (1.87383, 0.87), 
        (1.89223, 0.88), (1.91084, 0.89), (1.92964, 0.9), (1.94862, 0.91), (
        1.96775, 0.92), (1.98702, 0.93), (2.00638, 0.94), (2.02581, 0.95), (
        2.04527, 0.96), (2.06472, 0.97), (2.08411, 0.98), (2.1034, 0.99), (
        2.12253, 1.0)))
    createdModel.materials['VLGQ6590BAY'].Density(table=((1.2e-09, 
        ), ))
    
    # BONN-60SHA
    createdModel.Material(name='BONN-60SHA')
    createdModel.materials['BONN-60SHA'].Hyperelastic(
        materialType=const.ISOTROPIC, testData=const.OFF, type=const.MOONEY_RIVLIN, 
        volumetricResponse=const.VOLUMETRIC_DATA, table=((0.42981, 0.10745, 0.005), 
        ))
    createdModel.materials['BONN-60SHA'].Density(table=((1.0707E-09, 
        ), ))