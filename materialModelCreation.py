import numpy as np
import abaqusConstants as const
import os
import csv

class CalculateMaterialProperties:
    def __init__(self,Names) :
        self.symmtery_error = []
        self.homogenization_error = {}
        self.stiffness = {}
        parent1 = 'stiffnessMatrices'
        parent2 = Names.JOB_NAME
        self.parent_folder = os.path.join(parent1,parent2)
        try: 
            if not os.path.exists(self.parent_folder):
                os.makedirs(self.parent_folder)
        except OSError:
            if not os.path.isdir(self.parent_folder):
                raise
    def calculateMatrixDiffNorm(self,Matrix1, Matrix2):
        """ Return the relative norm of the difference between two matrices"""
        return np.linalg.norm(np.subtract(Matrix1, Matrix2)) / np.linalg.norm(Matrix2)
    def getSymmetricMatrix(self,fullMatrix):
        """
        Return a symmetrized 2D numpy array (with corresponding off-diagonal terms averaged)
        """
        if not isinstance(fullMatrix, np.ndarray):
            fullMatrix = np.array(fullMatrix)
        symmMatrix = np.array(object=fullMatrix)
        shape = fullMatrix.shape
        if not (len(shape) == 2 and shape[0] == shape[1]):
            raise TypeError('input matrix is not 2D and square, cannot symmetrize')
        for i in range(shape[0]):
            for j in range(i):
                avg = (fullMatrix[(i, j)] + fullMatrix[(j, i)]) / 2.0
                symmMatrix[(i, j)] = symmMatrix[(j, i)] = avg

        return symmMatrix
    def addIsoElasticityDescription(self,material, cMatrix):
        """
        Approximate the stiffness matrix as an isotropic material description
        and add it to the material.
        """
        cSymm = self.getSymmetricMatrix(cMatrix)
        symmetricRelativeErrorNorm = self.calculateMatrixDiffNorm(cSymm, cMatrix)
        symmErrorLine = 'Symmetry Error norm((C+C^T)/2.0 - C)/norm(C) = %12.5g%%' % (symmetricRelativeErrorNorm * 100)
        print(symmErrorLine)
        self.symmtery_error = symmetricRelativeErrorNorm
        
        mu1 = (cSymm[3][3] + cSymm[4][4] + cSymm[5][5]) / 3.0
        l1 = (cSymm[0][1] + cSymm[0][2] + cSymm[1][2]) / 3.0
        l2 = (cSymm[0][0] + cSymm[1][1] + cSymm[2][2]) / 3.0 - 2.0 * mu1
        mu2 = (cSymm[0][0] + cSymm[1][1] + cSymm[2][2]) / 6.0 - 0.5 * l1
        L = (l1 + l2) / 2.0
        mu = 1.0 / 2.0 * (mu1 + mu2)
        E = mu * (3.0 * L + 2.0 * mu) / (L + mu)
        poissonRatio = L / (2.0 * (L + mu))
        material.Elastic(table=((E, poissonRatio),))
        coeff = E / (1 + poissonRatio) / (1 - 2 * poissonRatio)
        C_iso = np.zeros((6, 6))
        C_iso[0:3, 0:3] = coeff * poissonRatio
        C_iso[0][0] = coeff * (1 - poissonRatio)
        C_iso[1][1] = coeff * (1 - poissonRatio)
        C_iso[2][2] = coeff * (1 - poissonRatio)
        C_iso[3][3] = coeff * (1 - 2 * poissonRatio) / 2.0
        C_iso[4][4] = coeff * (1 - 2 * poissonRatio) / 2.0
        C_iso[5][5] = coeff * (1 - 2 * poissonRatio) / 2.0
        IsoRelativeErrorNorm = self.calculateMatrixDiffNorm(C_iso, cSymm)
        self.homogenization_error['Isotropy'] = IsoRelativeErrorNorm
        self.stiffness['isotropy'] = C_iso
        file_name = os.path.join(self.parent_folder,'isotropic-StiffnessMatrix.txt')
        np.savetxt(file_name, C_iso, delimiter=',', fmt='%1.4e')
        if IsoRelativeErrorNorm > 0.1:
            print('Isotropy assumed in ' + material.name + ' may not be a good assumption')
            print('The isotropic stiffness matrix differs by more than 10% from the')
            print('full anisotropic stiffness matrix')
        IsoErrorLine = 'Isotropy Error norm(C_iso-C_symm)/norm(C_symm) = %12.5g%%' % (IsoRelativeErrorNorm * 100)
        print(IsoErrorLine)
    def addEngineeringConstantElasticityDescription(self,material, cMatrix):
        """
        Approximate the stiffness matrix as an orthotropic material description
        and add it to the material using an engineering constnats elasticity
        definition
        """
        cSymm = self.getSymmetricMatrix(cMatrix)
        symmetricRelativeErrorNorm = self.calculateMatrixDiffNorm(cSymm, cMatrix)
        symmErrorLine = 'Symmetry Error norm((C+C^T)/2.0 - C)/norm(C) = %12.5g%%' % (symmetricRelativeErrorNorm * 100)
        print(symmErrorLine)
        self.symmtery_error = symmErrorLine
        
        C_ortho = np.zeros((6, 6))
        C_ortho[0:3, 0:3] = cSymm[0:3, 0:3]
        C_ortho[3][3] = cSymm[3][3]
        C_ortho[4][4] = cSymm[4][4]
        C_ortho[5][5] = cSymm[5][5]
        S = np.linalg.inv(cSymm)
        E1 = 1.0 / S[0][0]
        E2 = 1.0 / S[1][1]
        E3 = 1.0 / S[2][2]
        material.Elastic(type=const.ENGINEERING_CONSTANTS, table=(
        (
        E1, E2, E3,
        -S[0][1] * E1, -S[0][2] * E1, -S[1][2] * E2,
        1.0 / S[3][3], 1.0 / S[4][4], 1.0 / S[5][5]),))
        OrthoRelativeErrorNorm = self.calculateMatrixDiffNorm(C_ortho, cSymm)
        OrthoErrorLine = 'Orthotropy Error norm(C_ortho-C_symm)/norm(C_symm) = %12.5g%%' % (OrthoRelativeErrorNorm * 100)
        print(OrthoErrorLine)
        self.homogenization_error['Otrhotropy_eng_constants'] = OrthoRelativeErrorNorm
        
        if OrthoRelativeErrorNorm > 0.1:
            print('Orthotropy assumed in ' + material.name + ' may not be a good assumption')
            print('The orthotropic stiffness matrix differs by more than 10% from the')
            print('full anisotropic stiffness matrix')


    def addOrthotropicElasticityDescription(self,material, cMatrix):
        """
        Approximate the stiffness matrix as an orthotropic material description
        and add it to the material using an orthotropic elasticity
        definition
        """
        cSymm = self.getSymmetricMatrix(cMatrix)
        symmetricRelativeErrorNorm = self.calculateMatrixDiffNorm(cSymm, cMatrix)
        symmErrorLine = 'Symmetry Error norm((C+C^T)/2.0 - C)/norm(C) = %12.5g%%' % (symmetricRelativeErrorNorm * 100)
        print(symmErrorLine)
        self.symmtery_error = symmErrorLine
        
        C_ortho = np.zeros((6, 6))
        C_ortho[0:3, 0:3] = cSymm[0:3, 0:3]
        C_ortho[3][3] = cSymm[3][3]
        C_ortho[4][4] = cSymm[4][4]
        C_ortho[5][5] = cSymm[5][5]
        material.Elastic(type=const.ORTHOTROPIC, table=(
        (
        cSymm[0][0], cSymm[0][1], cSymm[1][1],
        cSymm[0][2], cSymm[1][2], cSymm[2][2],
        cSymm[3][3], cSymm[4][4], cSymm[5][5]),))
        OrthoRelativeErrorNorm = self.calculateMatrixDiffNorm(C_ortho, cMatrix)
        OrthoErrorLine = 'Orthotropy Error norm(C_ortho-C_symm)/norm(C_symm) = %12.5g%%' % (OrthoRelativeErrorNorm * 100)
        print(OrthoErrorLine)
        file_name = os.path.join(self.parent_folder,'orthotropic-StiffnessMatrix.txt')
        np.savetxt(file_name, C_ortho, delimiter=',', fmt='%1.4e')
        self.homogenization_error['Otrhotropy'] = OrthoRelativeErrorNorm
        self.stiffness['orthotropy'] = C_ortho
        
        if OrthoRelativeErrorNorm > 0.1:
            print('Orthotropy assumed in ' + material.name + ' may not be a good assumption')
            print('The orthotropic stiffness matrix differs by more than 10% from the')
            print('full anisotropic stiffness matrix')
    def addTransverseElasticityDescription(self,material, cMatrix, symmertryAxis = 'z'):
        """
        Write the stiffness matrix to the material as an Transverse definition
        """
        cSymm = self.getSymmetricMatrix(cMatrix)
        symmetricRelativeErrorNorm = self.calculateMatrixDiffNorm(cSymm, cMatrix)
        symmErrorLine = 'Symmetry Error norm((C+C^T)/2.0 - C)/norm(C) = %12.5g%%' % (symmetricRelativeErrorNorm * 100)
        print(symmErrorLine)
        self.symmtery_error = symmErrorLine
        C_t = np.zeros((6, 6))
        C = {}
        # C_mat = np.zeros((6,6))
        for i, row_num in enumerate(range(1,7)):
            indices = [str(row_num)+str(col_num) for col_num in range(1,7)]
            for idx, ind in enumerate(indices):
                C[ind] = float(cSymm[i,idx])
        # TODO: print these engineering constants
        EL = C['11']-2*C['12']*C['21']/(C['22']+C['23'])  # Longitudinal E1 modulus
        nuL = C['12']/(C['22']+C['23'])  # 12 Poisson coefficient
        ET = (C['11']*(C['22']+C['23'])-2*C['12']*C['12']) * \
            (C['22']-C['23'])/(C['11']*C['22']-C['12']*C['21'])  # E2
        nuT = (C['11']*C['23']-C['12']*C['21']) / \
            (C['11']*C['22']-C['12']*C['21'])  # 23 Poisson coefficient
        GT = (C['22']-C['23'])/2.0  # or GT=ET/2/(1+nuT) # 23 Shear stiffness
        eng_constants = {'E_L':EL, 'poisson_L':nuL,
                         'E_T':ET,'poisson_T':nuT,
                         'G_T':GT}
        print("Calculated engineering constants for the transverse isotropy materail:")
        for k,v in eng_constants.items():
            msg = k + ' = ' + str(v)
            print(msg)
        C_t[0,0] = C_t[1,1]=C['11']
        C_t[0,1] = C_t[1,0]=C['12']
        for i,j in zip([0,1,2,2],[2,2,0,1]):
            C_t[i,j] = C['13']
        coef1 = (C['11']-C['12'])/2.0
        if symmertryAxis.lower() == 'z':
            coef2 = C['44']
            C_t[3,3] = coef1
            C_t[5,5] = C_t[4,4] = coef2
        if symmertryAxis.lower() == 'y':
            coef2 = C['55']
            C_t[4,4] = coef1
            C_t[5,5] = C_t[3,3] = coef2
        if symmertryAxis.lower() == 'x':
            coef2 = C['66']
            C_t[5,5] = coef1
            C_t[3,3] = C_t[4,4] = coef2
        TransverseRelativeErrorNorm = self.calculateMatrixDiffNorm(C_t, cMatrix)
        TransverseErrorLine = 'Transverse Error norm(C_t-C_symm)/norm(C_symm) = %12.5g%%' % (TransverseRelativeErrorNorm * 100)
        print(TransverseErrorLine)
        self.homogenization_error['Transverse_'+symmertryAxis] = TransverseRelativeErrorNorm
        self.stiffness['Transverse_'+symmertryAxis] = C_t
        file_name = os.path.join(self.parent_folder,'transverse-%s-StiffnessMatrix.txt'% (symmertryAxis))
        np.savetxt(file_name , C_t, delimiter=',', fmt='%1.4e')
        
        if TransverseRelativeErrorNorm > 0.1:
            print('Transerve assumed in ' + material.name + ' may not be a good assumption')
            print('The transverse stiffness matrix differs by more than 10% from the')
            print('full anisotropic stiffness matrix')



    
    def addMonotropicElasticityDescription(self,material, cMatrix):
        """
        Write the stiffness matrix to the material as an monotropic elasticity definition
        """
        cSymm = self.getSymmetricMatrix(cMatrix)
        symmetricRelativeErrorNorm = self.calculateMatrixDiffNorm(cSymm, cMatrix)
        symmErrorLine = 'Symmetry Error norm((C+C^T)/2.0 - C)/norm(C) = %12.5g%%' % (symmetricRelativeErrorNorm * 100)
        print(symmErrorLine)
        self.symmtery_error = symmErrorLine
        C_mono = np.zeros((6, 6))
        C_mono[0:4, 0:4] = cSymm[0:4, 0:4]
        C_mono[4][4] = cSymm[4][4]
        C_mono[5][5] = cSymm[5][5]
        C_mono[5][4] = C_mono[4][5] = cSymm[4][5]

        material.Elastic(type=const.ANISOTROPIC, table=([ C_mono[i][j] for i in range(6) for j in range(i + 1) ],))

        monotropicRelativeErrorNorm = self.calculateMatrixDiffNorm(C_mono, cMatrix)
        MonoErrorLine = 'Monotropic Error norm(C_ortho-C_symm)/norm(C_symm) = %12.5g%%' % (monotropicRelativeErrorNorm * 100)
        print(MonoErrorLine)
        file_name = os.path.join(self.parent_folder,'monotropic-StiffnessMatrix.txt')
        np.savetxt(file_name, C_mono, delimiter=',', fmt='%1.4e')
        self.homogenization_error['Monotropy'] = monotropicRelativeErrorNorm
        self.stiffness['monotropy'] = C_mono
        if monotropicRelativeErrorNorm > 0.1:
            print('Monotropy assumed in ' + material.name + ' may not be a good assumption')
            print('The orthotropic stiffness matrix differs by more than 10% from the')
            print('full anisotropic stiffness matrix')

    
    def addAnisotropicElasticityDescription(self,material, cMatrix):
        """
        Write the stiffness matrix to the material as an anisotropic elasticity definition
        """
        cSymm = self.getSymmetricMatrix(cMatrix)
        symmetricRelativeErrorNorm = self.calculateMatrixDiffNorm(cSymm, cMatrix)
        symmErrorLine = 'Symmetry Error norm((C+C^T)/2.0 - C)/norm(C) = %12.5g%%' % (symmetricRelativeErrorNorm * 100)
        print(symmErrorLine)
        self.symmtery_error = symmErrorLine
        file_name = os.path.join(self.parent_folder,'anisotropic-StiffnessMatrix.txt')
        np.savetxt(file_name, cSymm, delimiter=',', fmt='%1.4e')
        self.stiffness['anisotropy'] = cSymm
        
        material.Elastic(type=const.ANISOTROPIC, table=([ cSymm[i][j] for i in range(6) for j in range(i + 1) ],))

def generateMaterialsMechanical(model, nameSuffix, Names, cMatrix=None, aniso=False, ortho=False, iso=False, engConst=False,mono=False, transverse=False, transverse_axes=['z']):
    """
    Create requested types of material objects in the specified CAE model
    
    Inputs:
    model -       model object in which materials are generated
    nameSuffix -  string that will idenfity which step this homogenization was performed for
    Names      -  class which keeps the global information
    cMatrix -     stiffness matrix.  If provided, elasticity definitions will be created
    aniso/ortho/iso/engConst - bools specifying what types of definitions to create
    """
    calMatProp = CalculateMaterialProperties(Names)
    checkCMatrix = not any(None in sub for sub in cMatrix)
    if checkCMatrix:
        cSymm = calMatProp.getSymmetricMatrix(cMatrix)
        symmetricRelativeErrorNorm = calMatProp.calculateMatrixDiffNorm(cSymm, cMatrix)
        if symmetricRelativeErrorNorm > 0.1:
            print('The stiffness calculated for materials with the suffix ' + nameSuffix + ' is based')
            print('on a stiffness matrix that differs from its symmetrized version by %f6.3%%' % (symmetricRelativeErrorNorm * 100))
    
    if iso:
        isoMat = model.Material(name='Homog_Iso_' + nameSuffix)
        if checkCMatrix:
            calMatProp.addIsoElasticityDescription(isoMat, cMatrix)
        
        
    if engConst:
        engConstMat = model.Material(name='Homog_EngConst_' + nameSuffix)
        if checkCMatrix:
            calMatProp.addEngineeringConstantElasticityDescription(engConstMat, cMatrix)
        
    if ortho:
        orthoMat = model.Material(name='Homog_Ortho_' + nameSuffix)
        if checkCMatrix:
            calMatProp.addOrthotropicElasticityDescription(orthoMat, cMatrix)
        
    if aniso:
        anisoMat = model.Material(name='Homog_Aniso_' + nameSuffix)
        if checkCMatrix:
            calMatProp.addAnisotropicElasticityDescription(anisoMat, cMatrix)
    if mono:
        monoMat = model.Material(name='Homog_Mono_' + nameSuffix)
        if checkCMatrix:
            calMatProp.addMonotropicElasticityDescription(monoMat, cMatrix)
    if transverse:
        if checkCMatrix:
            for transverse_axis in transverse_axes:
                transverseMat = model.Material(name='Homog_transverse_' + transverse_axis + '_' + nameSuffix)
                calMatProp.addTransverseElasticityDescription(transverseMat, cMatrix,transverse_axis)
        
    return calMatProp

def generateAllHomogenizationScenerios(model, nameSuffix, Names, cMatrix=None):
    """
    Create all types of material objects in the specified CAE model
    
    Inputs:
    model -       model object in which materials are generated
    nameSuffix -  string that will idenfity which step this homogenization was performed for
    Names      -  class which keeps the global information
    cMatrix -     stiffness matrix.  If provided, elasticity definitions will be create
    """
    aniso=True; ortho=True; iso=True; engConst=False;mono=True; transverse=True; 
    transverse_axes=['x','y','z']
    matProp = generateMaterialsMechanical(model=model, nameSuffix=nameSuffix, 
                                                       Names=Names, 
                                                       cMatrix=cMatrix, 
                                                       aniso=aniso, 
                                                       ortho=ortho, 
                                                       iso=iso, 
                                                       engConst=engConst,
                                                       mono=mono, 
                                                       transverse=transverse, 
                                                       transverse_axes=transverse_axes)
    with open(os.path.join(matProp.parent_folder,'%s_%s_homogenzation_erros.csv'%(Names.JOB_NAME,nameSuffix)),'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        for k,v in matProp.homogenization_error.items():
            csv_writer.writerow([k,v])
        csv_writer.writerow(['Symmetry error',matProp.symmtery_error])
    return matProp