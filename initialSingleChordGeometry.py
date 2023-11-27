
from simulationSettingsMicro import Names
import abaqusConstants as const


def createInitialSketch(createdModel,a1,a2,a3):
    createdModel.ConstrainedSketch(name='square_cross_section', sheetSize=10.0)
    createdModel.sketches['square_cross_section'].rectangle(point1=(-a1, -a2), 
        point2=(a1, a2))
    createdModel.Part(dimensionality=const.THREE_D, name=Names.PART_NAME, type=
        const.DEFORMABLE_BODY)
    createdModel.parts[Names.PART_NAME].BaseSolidExtrude(depth=2*a3, sketch=
        createdModel.sketches['square_cross_section'])
    del createdModel.sketches['square_cross_section']


# partition of the fiber
def create_partition(createdModel,a1,a2,a3,r_f):
    p = createdModel.parts[Names.PART_NAME]
    f, e, d1 = p.faces, p.edges, p.datums
    t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(a1, 0.0, 
        a3)), sketchUpEdge=e.findAt(coordinates=(a1, 2*a1, 
        0.0)), sketchPlaneSide=const.SIDE1, origin=(a1, 0.0, a3))
    s = createdModel.ConstrainedSketch(name='__profile__', sheetSize=2.43, 
        gridSpacing=0.06, transform=t)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=const.SUPERIMPOSE)
    p.projectReferencesOntoSketch(sketch=s, filter=const.COPLANAR_EDGES)
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(r_f, 0.0))
    p = createdModel.parts[Names.PART_NAME]
    f = p.faces
    pickedFaces = f.findAt(((a1, 0.0, a3), ))
    e1, d2 = p.edges, p.datums
    p.PartitionFaceBySketch(sketchUpEdge=e1.findAt(coordinates=(a1, 2*a1, 
        0.0)), faces=pickedFaces, sketch=s)
    s.unsetPrimaryObject()
    del createdModel.sketches['__profile__']
    p = createdModel.parts[Names.PART_NAME]
    c = p.cells
    pickedCells = c.findAt(((-a1, 0.0, a3), ))
    e, d1 = p.edges, p.datums
    pickedEdges =(e.findAt(coordinates=(a1, r_f, a3)), )
    p.PartitionCellByExtrudeEdge(line=e.findAt(coordinates=(0.0, a2, 
        2*a3)), cells=pickedCells, edges=pickedEdges, sense=const.REVERSE)




def pipeLineInitialGeometry(createdModel,a1,a2,a3,r_f):
    createInitialSketch(createdModel=createdModel,a1=a1,a2=a2,a3=a3)
    create_partition(createdModel=createdModel,a1=a1,a2=a2,a3=a3,r_f=r_f)