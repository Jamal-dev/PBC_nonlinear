import math
import abaqusConstants as const
import numpy as np

def addingDepthPipeLine(depth,a1,a2,a3,r_f,s1):
    def circle_line_segment_intersection(pt1, pt2,circle_center, radius,  full_line=False, tangent_tol=1e-9):
        """ Find the points at which a circle intersects a line-segment.  This can happen at 0, 1, or 2 points.
        :param full_line: True to find intersections along full line - not just in the segment.  False will just return intersections within the segment.
        :param tangent_tol: Numerical tolerance at which we decide the intersections are close enough to consider it a tangent
        :return Sequence[Tuple[float, float]]: A list of length 0, 1, or 2, where each element is a point at which the circle intercepts a line segment.
        Note: We follow: http://mathworld.wolfram.com/Circle-LineIntersection.html
        """
        (p1x, p1y), (p2x, p2y), (cx, cy) = pt1, pt2, circle_center
        (x1, y1), (x2, y2) = (p1x - cx, p1y - cy), (p2x - cx, p2y - cy)
        dx, dy = (x2 - x1), (y2 - y1)
        dr = (dx ** 2 + dy ** 2)**.5
        big_d = x1 * y2 - x2 * y1
        discriminant = radius ** 2 * dr ** 2 - big_d ** 2
        if discriminant < 0:  # No intersection between circle and line
            return []
        else:  # There may be 0, 1, or 2 intersections with the segment
            intersections = [
                (cx + (big_d * dy + sign * (-1 if dy < 0 else 1) * dx * discriminant**.5) / dr ** 2,
                cy + (-big_d * dx + sign * abs(dy) * discriminant**.5) / dr ** 2)
                for sign in ((1, -1) if dy < 0 else (-1, 1))]  # This makes sure the order along the segment is correct
            if not full_line:  # If only considering the segment, filter out intersections that do not fall within the segment
                fraction_along_segment = [(xi - p1x) / dx if abs(dx) > abs(dy) else (yi - p1y) / dy for xi, yi in intersections]
                intersections = [pt for pt, frac in zip(intersections, fraction_along_segment) if 0 <= frac <= 1]
            if len(intersections) == 2 and abs(discriminant) <= tangent_tol:  # If line is tangent to circle, return just one point (as both intersections have same location)
                return [intersections[0]]
            else:
                return intersections
        
    def find_by_intersection(p1,p2,radius=r_f,z_c=a3,y_c=0.0):
        line_start = np.array([p1[0],p1[1]])
        line_end   = np.array([p2[0],p2[1]])
        circle_center = np.array([y_c,z_c])
        l = circle_line_segment_intersection(line_start, line_end, circle_center=circle_center, radius=radius)
        return (l[0][0],l[0][1])
    
    def mid_point(p1,p2):
        return tuple([(x1+x2)/2.0 for x1,x2 in zip(p1,p2)])
    
    def calculateAllVertices(a2=a2,a3=a3,r_f=r_f):
        r_i = r_f/2.0
        # the numbering is anticlockwise
        vertices_hex = []
        ang = math.cos(math.pi/4.0)
        vertices_hex.append((0.0,a3-r_i*ang))
        vertices_hex.append((r_i*ang,a3-r_i*ang))
        p1 = (r_i*ang,a3-r_i*ang)
        p2 = (r_i,a3)
        vertices_hex.append(mid_point(p1,p2))
        vertices_hex.append(p2)
        p1 = p2
        p2 = (r_i*ang,a3+r_i*ang)
        vertices_hex.append(mid_point(p1,p2))
        vertices_hex.append(p2)
        vertices_hex.append((0.0,a3+r_i*ang))
        p1 = (-r_i*ang,a3+r_i*ang)
        vertices_hex.append(p1)
        p2 = (-r_i,a3) 
        vertices_hex.append(mid_point(p1,p2))
        vertices_hex.append(p2)
        p1 = p2
        p2 = (-r_i*ang,a3-r_i*ang)
        vertices_hex.append(mid_point(p1,p2))
        vertices_hex.append(p2)
        # vertices of the outer most geometry
        vertices_outer = []
        vertices_outer.append((0.0,0.0))
        vertices_outer.append((a2/2.0,0.0))
        vertices_outer.append((a2,0.0))
        vertices_outer.append((a2,a3))
        vertices_outer.append((a2,2*a3))
        vertices_outer.append((a2/2.0,2*a3))
        vertices_outer.append((0.0,2*a3))
        vertices_outer.append((-a2/2.0,2*a3))
        vertices_outer.append((-a2,2*a3))
        vertices_outer.append((-a2,a3))
        vertices_outer.append((-a2,0.0))
        vertices_outer.append((-a2/2.0,0.0))
        # vertices at the circle
        vertices_circle = []
        vertices_circle.append((0.0,a3-r_f))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[1],p2=vertices_outer[1]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[2],p2=vertices_outer[2]))
        vertices_circle.append((r_f,a3))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[4],p2=vertices_outer[4]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[5],p2=vertices_outer[5]))
        vertices_circle.append((0.0,r_f + a3))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[7],p2=vertices_outer[7]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[8],p2=vertices_outer[8]))
        vertices_circle.append((-r_f, a3))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[10],p2=vertices_outer[10]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[11],p2=vertices_outer[11]))
        # vertex origin
        vertex_origin = (0.0,a3)
        vertices={}
        vertices['hex'] = vertices_hex
        vertices['outer'],vertices['circle'], = vertices_outer,vertices_circle
        vertices['origin'] = (vertex_origin[1],vertex_origin[0])
        keys = ['hex','outer','circle']
        for k in keys:
            v = vertices[k]
            temp = []
            for coord_y_z in v:
                temp.append((coord_y_z[1],coord_y_z[0]))
            vertices[k] = temp
        return vertices
    
    
    def addMidPoints(vertices):
        nextList = vertices[1:] + [vertices[0]]
        newList = []
        for index,(curPoint, nextPoint) in enumerate(zip(vertices,nextList)):
            midPoint = mid_point(curPoint,nextPoint)
            newList.append(curPoint)
            newList.append(midPoint)
        return newList
    
    def addDepth(d,vertices):
        points = vertices
        for _ in range(d):
            points = addMidPoints(points)
        return points
    
    def set_difference(new_list, previous_list):
        def is_close(a, b, rel_tol=1e-9, abs_tol=0.0):
            return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
        def tuple_close(t1, t2):
            if len(t1) != len(t2):
                return False
            for a, b in zip(t1, t2):
                if not is_close(a, b):
                    return False
            return True
        result = []
        for item in new_list:
            if not any(tuple_close(item, prev_item) for prev_item in previous_list):
                result.append(item)
        return result
    
    def popPrevVertices(newVertices,prevVertices):
        newDic = {}
        keys = ['hex','outer','circle']
        for k in keys:
            # newDic[k] = list(set(newVertices[k]) - set(prevVertices[k]))
            newDic[k] = set_difference(newVertices[k],prevVertices[k]) 
        newDic['origin'] = prevVertices['origin']
        return newDic
    
    def pipeLineAddingDepth(depth, prevVertices_dic):
        keys = ['hex','outer']
        newVertices_dic = {}
        for k in keys:
            v = prevVertices_dic[k]
            newVertices_dic[k] = addDepth(depth,v)
        x_c, y_c = prevVertices_dic['origin']
        newVertices_dic['origin'] = prevVertices_dic['origin']
        newVertices_dic['circle'] = []
        for p1,p2 in zip(newVertices_dic[keys[0]],newVertices_dic[keys[1]]):
            newVertices_dic['circle'].append(find_by_intersection(p1=p1,p2=p2,z_c = y_c, y_c = x_c))
        cleanedVertices_dic = popPrevVertices(newVertices=newVertices_dic,prevVertices=prevVertices_dic)
        return cleanedVertices_dic, newVertices_dic
    
    def reverseCoordinates(vertices,appendingCoord=a1):
        for k,v in vertices.items():
            if not 'origin' in k:
                temp = []
                for coord_z_y in v:
                    temp.append((appendingCoord,coord_z_y[1],coord_z_y[0]))
                vertices[k] = temp
            else:
                vertices[k] = (appendingCoord,v[1],v[0])
        return vertices
    
    vertices = calculateAllVertices(a2=a2,a3=a3,r_f=r_f)
    cleanedVertices, newVertices = pipeLineAddingDepth(depth,vertices)
    newVertices3d = reverseCoordinates(newVertices,appendingCoord=a1)
    def shiftXcoord(c,shift=-a3):
        return (c[0]+shift,c[1])
    def drawLines(vertices,s1):
        keys = ['hex','circle','outer']
        for c1,c2,c3 in zip(vertices[keys[0]],vertices[keys[1]],vertices[keys[2]]):
            c1,c2,c3 = shiftXcoord(c1), shiftXcoord(c2), shiftXcoord(c3)
            s1.Line(point1=c1, point2=c2)
            s1.Line(point1=c2, point2=c3)
    drawLines(cleanedVertices,s1)
    return newVertices3d

def createSketchSingleChord(Model,r_f,a1,a2,a3,Names,depth):
    sheatSize = max(2*a2,2*a3) + max(2*a2,2*a3)/2.0 
    diag_sheat_size = sheatSize * math.cos(math.pi/4.0)   
    s1 = Model.ConstrainedSketch(name='__profile__', sheetSize=sheatSize)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=const.STANDALONE)
    eps = 1e-4
    r_i = r_f/2.0
    # creating construction line
    s1.ConstructionLine(point1=(0.0, 0.0), angle=0.0)
    s1.HorizontalConstraint(entity=g.findAt((sheatSize, 0.0)), addUndoState=False)
    s1.ConstructionLine(point1=(0.0, 0.0), angle=90.0)
    s1.VerticalConstraint(entity=g.findAt((0.0, sheatSize)), addUndoState=False)
    s1.ConstructionLine(point1=(0.0, 0.0), angle=45.0)
    s1.ConstructionLine(point1=(0.0, 0.0), angle=-45.0)
    s1.ConstructionCircleByCenterPerimeter(center=(0.0, 0.0), point1=(r_f, 0.0))
    s1.ConstructionCircleByCenterPerimeter(center=(0.0, 0.0), point1=(r_f/2.0, 0.0))
    # vertical line above line
    s1.Line(point1=(0.0, 0.0), point2=(0.0, a2))
    s1.VerticalConstraint(entity=g.findAt((0.0, (a2-eps))), addUndoState=False)
    # the geometry entry at (sheatSize,0 is the construction axis)
    s1.PerpendicularConstraint(entity1=g.findAt((sheatSize, 0.0)), entity2=g.findAt((0.0, 
        (a2-eps))), addUndoState=False)
    s1.CoincidentConstraint(entity1=v.findAt((0.0, 0.0)), entity2=g.findAt((sheatSize, 
        0.0)), addUndoState=False)
    # vertical line down line
    s1.Line(point1=(0.0, 0.0), point2=(0.0, -a2))
    s1.VerticalConstraint(entity=g.findAt((0.0, -(a2-eps))), addUndoState=False)
    s1.ParallelConstraint(entity1=g.findAt((0.0, (a2-eps))), entity2=g.findAt((0.0, 
        -(a2-eps))), addUndoState=False)
    # horizontal line
    s1.Line(point1=(0.0, 0.0), point2=(a3, 0.0))
    s1.HorizontalConstraint(entity=g.findAt(((a3-eps), 0.0)), addUndoState=False)
    s1.PerpendicularConstraint(entity1=g.findAt((0.0, (a2-eps))), entity2=g.findAt(
        ((a3-eps), 0.0)), addUndoState=False)
    # horizontal negative line
    s1.Line(point1=(0.0, 0.0), point2=(-a3, 0.0))
    s1.HorizontalConstraint(entity=g.findAt((-(a3-eps), 0.0)), addUndoState=False)
    s1.PerpendicularConstraint(entity1=g.findAt((0.0, (a2-eps))), entity2=g.findAt(
        (-(a3-eps), 0.0)), addUndoState=False)
    # upper left line
    s1.Line(point1=(-a3, 0.0), point2=(-a3, a2))
    s1.VerticalConstraint(entity=g.findAt((-a3, a2/2.0)), 
        addUndoState=False)
    s1.PerpendicularConstraint(entity1=g.findAt((-(a3-eps), 0.0)), 
        entity2=g.findAt((-a3, a2/2.0)), addUndoState=False)
    s1.Line(point1=(-a3, a2), point2=(0.0, 
        a2))
    s1.HorizontalConstraint(entity=g.findAt((-a3/2.0, a2)), 
        addUndoState=False)
    s1.PerpendicularConstraint(entity1=g.findAt((-a3, a2/2.0)), 
        entity2=g.findAt((-a3/2.0, a2)), addUndoState=False)

    s1.Line(point1=(0.0, -a2), point2=(-a3, 
            -a2))
    s1.HorizontalConstraint(entity=g.findAt((-(a3/2.0), -a2)), 
        addUndoState=False)
    s1.PerpendicularConstraint(entity1=g.findAt((0.0, -(a2-eps))), 
        entity2=g.findAt((-(a3/2.0), -a2)), addUndoState=False)

    s1.Line(point1=(-a3, -a2), point2=(-a3, 0.0))
    s1.VerticalConstraint(entity=g.findAt((-a3, -a2/2.0)), 
        addUndoState=False)
    # s1.PerpendicularConstraint(entity1=g.findAt((-(a3/2.0), -a2)), 
    #     entity2=g.findAt((-a3, -a2/2.0)), addUndoState=False)

    s1.Line(point1=(a3, 0.0), point2=(a3, a2))
    s1.VerticalConstraint(entity=g.findAt((a3, a2/2.0)), 
        addUndoState=False)
    # s1.PerpendicularConstraint(entity1=g.findAt((0.350606, 0.0)), entity2=g.findAt(
    #     (a3, a2/2.0)), addUndoState=False)

    s1.Line(point1=(a3, a2), point2=(0.0, a2))
    s1.HorizontalConstraint(entity=g.findAt(((a3/2.0), a2)), 
        addUndoState=False)
    # s1.PerpendicularConstraint(entity1=g.findAt((a3, a2/2.0)), 
    #     entity2=g.findAt(((a3/2.0), a2)), addUndoState=False)

    s1.Line(point1=(a3, 0.0), point2=(a3, -a2))
    # s1.VerticalConstraint(entity=g.findAt((a3, -a2/2.0)), 
    #     addUndoState=False)
    # s1.PerpendicularConstraint(entity1=g.findAt((0.350606, 0.0)), entity2=g.findAt(
    #     (a3, -a2/2.0)), addUndoState=False)

    s1.Line(point1=(a3, -a2), point2=(0.0, 
        -a2))
    s1.HorizontalConstraint(entity=g.findAt(((a3/2.0), -a2)), 
        addUndoState=False)
    # s1.PerpendicularConstraint(entity1=g.findAt((a3, -a2/2.0)), 
    #     entity2=g.findAt(((a3/2.0), -a2)), addUndoState=False)

    # outer geometry is finish here

    def c2p(angle):
        angle = float(math.pi/180.0 * float(angle))
        return (r_i*math.cos(angle), r_i*math.sin(angle))     

    s1.Line(point1=c2p(45), point2=c2p(-45))
    s1.VerticalConstraint(entity=g.findAt((c2p(45)[0], (c2p(45)[0]-eps))), 
        addUndoState=False)
    s1.CoincidentConstraint(entity1=v.findAt((c2p(45)[0], -c2p(45)[0])), 
        entity2=g.findAt((diag_sheat_size, -diag_sheat_size)), addUndoState=False)

    s1.Line(point1=(-c2p(45)[0], -c2p(45)[0]), point2=(
        -c2p(45)[0], c2p(45)[0]))
    s1.VerticalConstraint(entity=g.findAt((-c2p(45)[0], -(c2p(45)[0]-eps))), 
        addUndoState=False)

    s1.Line(point1=(-c2p(45)[0], c2p(45)[0]), point2=(0.0, r_i))
    s1.CoincidentConstraint(entity1=v.findAt((0.0, r_i)), entity2=g.findAt((0.0, 
        sheatSize)), addUndoState=False)
    s1.Line(point1=(0.0, r_i), point2=(c2p(45)[0], c2p(45)[0]))

    s1.Line(point1=(-c2p(45)[0], -c2p(45)[0]), point2=(0.0, -r_i))
    s1.CoincidentConstraint(entity1=v.findAt((0.0, -r_i)), entity2=g.findAt((0.0, 
        sheatSize)), addUndoState=False)

    s1.Line(point1=(0.0, -r_i), point2=(c2p(45)[0], -c2p(45)[0]))

    def mid_point(a1,a2):
        p1 = c2p(a1)
        p2 = c2p(a2)
        return tuple([(x1+x2)/2.0 for x1,x2 in zip(p1,p2)])

    def pointOnLine(a1,a2):
        p1 = c2p(a1)
        p2 = c2p(a2)
        perturbPoint_x = p1[0] + 1e-2
        m = (p2[1]-p1[1])/(p2[0]-p1[0])
        b = p2[1] - m * p2[0]
        perturbPoint_y =  m * perturbPoint_x + b
        return (perturbPoint_x,perturbPoint_y)

    s1.Line(point1=mid_point(135,90), point2=(-a3, 
        a2))
    s1.CoincidentConstraint(entity1=v.findAt(mid_point(135,90)), 
        entity2=g.findAt(pointOnLine(135,90)), addUndoState=False)
    s1.EqualDistanceConstraint(entity1=v.findAt((-c2p(45)[0], c2p(45)[0])), 
        entity2=v.findAt((0.0, r_i)), midpoint=v.findAt(mid_point(135,90)), 
        addUndoState=False)

    s1.Line(point1=mid_point(45,90), point2=(a3, 
        a2))
    s1.CoincidentConstraint(entity1=v.findAt(mid_point(45,90)), 
        entity2=g.findAt(pointOnLine(90,45)), addUndoState=False)
    s1.EqualDistanceConstraint(entity1=v.findAt((0.0, r_i)), entity2=v.findAt((
        c2p(45)[0], c2p(45)[0])), midpoint=v.findAt(mid_point(45,90)), 
        addUndoState=False)

    s1.Line(point1=(-a3, a2/2.0), point2=(-c2p(45)[0], 
        c2p(45)[0]))
    s1.CoincidentConstraint(entity1=v.findAt((-a3, a2/2.0)), 
        entity2=g.findAt((-a3, (a2/2.0-eps))), addUndoState=False)
    s1.EqualDistanceConstraint(entity1=v.findAt((-a3, 0.0)), 
        entity2=v.findAt((-a3, a2)), midpoint=v.findAt((-a3, 
        a2/2.0)), addUndoState=False)

    s1.Line(point1=(c2p(45)[0], c2p(45)[0]), point2=(a3, 
        a2/2.0))
    s1.CoincidentConstraint(entity1=v.findAt((a3, a2/2.0)), 
        entity2=g.findAt((a3, (a2/2.0-eps))), addUndoState=False)
    s1.EqualDistanceConstraint(entity1=v.findAt((a3, 0.0)), entity2=v.findAt(
        (a3, a2)), midpoint=v.findAt((a3, a2/2.0)), 
        addUndoState=False)

    #################
    s1.Line(point1=mid_point(-90,-135), point2=(-a3, 
        -a2))
    s1.CoincidentConstraint(entity1=v.findAt(mid_point(-90,-135)), 
        entity2=g.findAt(pointOnLine(-135,-90)), addUndoState=False)

    s1.Line(point1=mid_point(-45,-90), point2=(a3, 
        -a2))
    s1.CoincidentConstraint(entity1=v.findAt(mid_point(-45,-90)), 
        entity2=g.findAt(pointOnLine(-90,-45)), addUndoState=False)

    s1.Line(point1=(-a3, -a2/2.0), point2=(-c2p(45)[0], 
        -c2p(45)[0]))
    s1.CoincidentConstraint(entity1=v.findAt((-a3, -a2/2.0)), 
        entity2=g.findAt((-a3, (-a2/2.0-eps))), addUndoState=False)

    s1.Line(point1=(c2p(45)[0], -c2p(45)[0]), point2=(a3, 
        -a2/2.0))
    s1.CoincidentConstraint(entity1=v.findAt((a3, -a2/2.0)), 
        entity2=g.findAt((a3, (-a2/2.0-eps))), addUndoState=False)
    Model.sketches.changeKey(fromName='__profile__', 
        toName=Names.PARTITION_FACE_SKETCH_NAME)
    updatedVertices = addingDepthPipeLine(depth=depth,a1=a1,a2=a2,a3=a3,r_f=r_f,s1=s1)
    s1.unsetPrimaryObject()
    return updatedVertices

def importSketchAndPartition(Model,a1,a2,a3,Names):
    a = Model.rootAssembly
    f1 = a.instances[Names.INSTANCE_NAME].faces
    e1 = a.instances[Names.INSTANCE_NAME].edges
    eps = 1e-4
    t = a.MakeSketchTransform(sketchPlane=f1.findAt(coordinates=(a1, 0.0, 2*a3-eps)), 
                              sketchUpEdge=e1.findAt(coordinates=(a1, 
        a2/2.0, 0.0)), sketchPlaneSide=const.SIDE1, origin=(a1, 0.0, 
        a3))
    s1 = Model.ConstrainedSketch(name='__profile__', 
        sheetSize=2.41, gridSpacing=0.06, transform=t)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=const.SUPERIMPOSE)
    a = Model.rootAssembly
    a.projectReferencesOntoSketch(sketch=s1, filter=const.COPLANAR_EDGES)
    s1.retrieveSketch(sketch=Model.sketches[Names.PARTITION_FACE_SKETCH_NAME])
    a = Model.rootAssembly
    f1 = a.instances[Names.INSTANCE_NAME].faces
    pickedFaces = f1.findAt(((a1, 0.0, 2*a3-eps), ), ((a1, 0.0, 
        a3), ))
    e11 = a.instances[Names.INSTANCE_NAME].edges
    a.PartitionFaceBySketch(sketchUpEdge=e11.findAt(coordinates=(a1, 
        a2/2.0, 0.0)), faces=pickedFaces, sketch=s1)
    s1.unsetPrimaryObject()
    del Model.sketches['__profile__']


def extrudePartionSketch(Model,r_f,a1,a2,a3,Names):
    a = Model.rootAssembly
    r_i = r_f/2.0
    e1 = a.instances[Names.INSTANCE_NAME].edges
    def mid_point(p1,p2,a1=a1,edge=True,e1 = e1):
        if len(p1) == 2:
            p1 = (a1,p1[1],p1[0])
        if len(p2) == 2:
            p2 = (a1,p2[1],p2[0])
        if edge:
            return e1.findAt((tuple([(x1+x2)/2.0 for x1,x2 in zip(p1,p2)]),))
        else:
            return tuple([(x1+x2)/2.0 for x1,x2 in zip(p1,p2)])


    def circle_line_segment_intersection(pt1, pt2,circle_center, radius,  full_line=False, tangent_tol=1e-9):
        """ Find the points at which a circle intersects a line-segment.  This can happen at 0, 1, or 2 points.
        :param full_line: True to find intersections along full line - not just in the segment.  False will just return intersections within the segment.
        :param tangent_tol: Numerical tolerance at which we decide the intersections are close enough to consider it a tangent
        :return Sequence[Tuple[float, float]]: A list of length 0, 1, or 2, where each element is a point at which the circle intercepts a line segment.
        Note: We follow: http://mathworld.wolfram.com/Circle-LineIntersection.html
        """
        (p1x, p1y), (p2x, p2y), (cx, cy) = pt1, pt2, circle_center
        (x1, y1), (x2, y2) = (p1x - cx, p1y - cy), (p2x - cx, p2y - cy)
        dx, dy = (x2 - x1), (y2 - y1)
        dr = (dx ** 2 + dy ** 2)**.5
        big_d = x1 * y2 - x2 * y1
        discriminant = radius ** 2 * dr ** 2 - big_d ** 2
        if discriminant < 0:  # No intersection between circle and line
            return []
        else:  # There may be 0, 1, or 2 intersections with the segment
            intersections = [
                (cx + (big_d * dy + sign * (-1 if dy < 0 else 1) * dx * discriminant**.5) / dr ** 2,
                cy + (-big_d * dx + sign * abs(dy) * discriminant**.5) / dr ** 2)
                for sign in ((1, -1) if dy < 0 else (-1, 1))]  # This makes sure the order along the segment is correct
            if not full_line:  # If only considering the segment, filter out intersections that do not fall within the segment
                fraction_along_segment = [(xi - p1x) / dx if abs(dx) > abs(dy) else (yi - p1y) / dy for xi, yi in intersections]
                intersections = [pt for pt, frac in zip(intersections, fraction_along_segment) if 0 <= frac <= 1]
            if len(intersections) == 2 and abs(discriminant) <= tangent_tol:  # If line is tangent to circle, return just one point (as both intersections have same location)
                return [intersections[0]]
            else:
                return intersections

    def find_by_intersection(p1,p2,radius=r_f,z_c=a3,y_c=0.0):
        line_start = np.array([p1[1],p1[2]])
        line_end = np.array([p2[1],p2[2]])
        circle_center = np.array([y_c,z_c])
        l = circle_line_segment_intersection(line_start, line_end, circle_center=circle_center, radius=radius)
        return (p1[0],l[0][0],l[0][1])

    def calculateAllVertices(a1=a1,a2=a2,a3=a3,r_f=r_f):
        # the numbering is anticlockwise
        vertices_hex = []
        ang = math.cos(math.pi/4.0)
        vertices_hex.append((a1,0.0,a3-r_i*ang))
        vertices_hex.append((a1,r_i*ang,a3-r_i*ang))
        p1 = (a1,r_i*ang,a3-r_i*ang)
        p2 = (a1,r_i,a3)
        vertices_hex.append(mid_point(p1,p2,edge=False))
        vertices_hex.append(p2)
        p1 = p2
        p2 = (a1,r_i*ang,a3+r_i*ang)
        vertices_hex.append(mid_point(p1,p2,edge=False))
        vertices_hex.append(p2)
        vertices_hex.append((a1,0.0,a3+r_i*ang))
        p1 = (a1,-r_i*ang,a3+r_i*ang)
        vertices_hex.append(p1)
        p2 = (a1,-r_i,a3) 
        vertices_hex.append(mid_point(p1,p2,edge=False))
        vertices_hex.append(p2)
        p1 = p2
        p2 = (a1,-r_i*ang,a3-r_i*ang)
        vertices_hex.append(mid_point(p1,p2,edge=False))
        vertices_hex.append(p2)
        # vertices of the outer most geometry
        vertices_outer = []
        vertices_outer.append((a1,0.0,0.0))
        vertices_outer.append((a1,a2/2.0,0.0))
        vertices_outer.append((a1,a2,0.0))
        vertices_outer.append((a1,a2,a3))
        vertices_outer.append((a1,a2,2*a3))
        vertices_outer.append((a1,a2/2.0,2*a3))
        vertices_outer.append((a1,0.0,2*a3))
        vertices_outer.append((a1,-a2/2.0,2*a3))
        vertices_outer.append((a1,-a2,2*a3))
        vertices_outer.append((a1,-a2,a3))
        vertices_outer.append((a1,-a2,0.0))
        vertices_outer.append((a1,-a2/2.0,0.0))
        # vertices at the circle
        vertices_circle = []
        vertices_circle.append((a1,0.0,a3-r_f))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[1],p2=vertices_outer[1]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[2],p2=vertices_outer[2]))
        vertices_circle.append((a1,r_f,a3))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[4],p2=vertices_outer[4]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[5],p2=vertices_outer[5]))
        vertices_circle.append((a1,0.0,r_f + a3))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[7],p2=vertices_outer[7]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[8],p2=vertices_outer[8]))
        vertices_circle.append((a1,-r_f, a3))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[10],p2=vertices_outer[10]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[11],p2=vertices_outer[11]))
        # vertex origin
        vertex_origin = (a1,0.0,a3)
        return vertices_hex,vertices_outer,vertices_circle,vertex_origin



    def arc_midpoint(start, end, center):
        if type(start) == tuple:
            start = np.array(list(start))
        if type(end) == tuple:
            end = np.array(list(end))
        if type(center) == tuple:
            center = np.array(list(center))
        # Calculate the vector from start to end point
        AB = end - start
        # Calculate the normal vector
        n = np.cross(AB, center - start)
        # Calculate the vector perpendicular to both AB and n
        v = np.cross(AB, n)
        # Normalize v to have length equal to the radius of the circular arc
        v = v / np.linalg.norm(v) * np.linalg.norm(center - start)
        # Calculate the midpoint
        midpoint = center + v
        return tuple(midpoint)

    def find_edges_by_vertices(list1,list2=None,midPointByArc=False,vertex_origin=vertex_origin,e1 = e1,a1=a1):
        if list2 is None:
            list2 = list1[1:]
            list2.append(list1[0])
        edges_list = []
        for p1,p2 in zip(list1,list2):
            if midPointByArc:
                midpoint = arc_midpoint(start=p2,end=p1,center=vertex_origin)
            else:
                midpoint = mid_point(p1,p2,a1=a1,edge=False,e1 = e1)
            edges_list.append(e1.findAt((midpoint,)))
        return edges_list

    def calculateAllEdges(vertices_hex,vertices_outer,vertices_circle,vertex_origin,e1=e1,a1=a1):
        edges_hex = find_edges_by_vertices(vertices_hex,midPointByArc=False,
                                        vertex_origin=vertex_origin,e1=e1,a1=a1)
        edges_circle = find_edges_by_vertices(vertices_circle,midPointByArc=True,
                                        vertex_origin=vertex_origin,e1=e1,a1=a1)
        edges_outer = find_edges_by_vertices(vertices_outer,midPointByArc=False,
                                        vertex_origin=vertex_origin,e1=e1,a1=a1)
        edges_staight_inner = find_edges_by_vertices(list1=vertices_hex,list2=vertices_circle,midPointByArc=False,
                                        vertex_origin=vertex_origin,e1=e1,a1=a1)
        edges_staight_outer = find_edges_by_vertices(list1=vertices_circle,list2=vertices_outer,midPointByArc=False,
                                        vertex_origin=vertex_origin,e1=e1,a1=a1)
        list1 = [vertex_origin for _ in range(4)]
        list2 = [vertices_hex[i] for i in [0,3,6,9]]
        edges_inside = find_edges_by_vertices(list1=list1,list2=list2,midPointByArc=False,
                                        vertex_origin=vertex_origin,e1=e1,a1=a1)
        return edges_hex,edges_circle,edges_outer,edges_staight_inner,edges_staight_outer,edges_inside



    def findCellEdges(edges1,edges2=None,edges3=None,edges4=None,num_edges_taken=3):
        special_scnerio = False
        if edges2 and edges4 is None:
            special_scnerio = True
        if not special_scnerio:
            if edges4 is None:
                edges4 = edges3[1:]
                edges4.append(edges3[0])
        if edges2 is None:
            edges2 = edges1[1:]
            edges2.append(edges2[0])
        cell_edges = []
        if not special_scnerio:
            for e1,e2,e3,e4 in zip(edges1,edges2,edges3,edges4):
                cell_edges.append(tuple([e1[0],e3[0],e2[0],e4[0]]))
        else:
            for i,(e1,e2) in enumerate(zip(edges1,edges2)):
                e3 = edges3[int(num_edges_taken*i):int(num_edges_taken*(i+1))]
                e3 = [edge[0] for edge in e3]
                e3.insert(0,e1[0])
                e3.extend([e2[0]])
                cell_edges.append(tuple(e3))
        return cell_edges

    def calculateAllCellEdges(edges_hex,edges_circle,edges_outer,edges_staight_inner,edges_staight_outer,edges_inside):
        cell_edges_level_1 = findCellEdges(edges_hex,edges_circle,edges_staight_inner)
        cell_edges_level_2 = findCellEdges(edges_circle,edges_outer,edges_staight_outer)
        cell_edges_level_0 = findCellEdges(edges1=edges_inside,edges2=None,edges3=edges_hex,edges4=None,num_edges_taken=3)
        return cell_edges_level_0,cell_edges_level_1,cell_edges_level_2

    vertices_hex,vertices_outer,vertices_circle,vertex_origin = calculateAllVertices(a1=a1,a2=a2,a3=a3,r_f=r_f)
    def calculateFreshPair(Names=Names,a1=a1,vertices_hex=vertices_hex,vertices_outer=vertices_outer,vertices_circle=vertices_circle,vertex_origin=vertex_origin):
        e1 = a.instances[Names.INSTANCE_NAME].edges
        edges_hex,edges_circle,edges_outer,edges_staight_inner,edges_staight_outer,edges_inside =calculateAllEdges(vertices_hex,vertices_outer,vertices_circle,vertex_origin,e1=e1,a1=a1)
        cell_edges_level_0,cell_edges_level_1,cell_edges_level_2=calculateAllCellEdges(edges_hex,edges_circle,edges_outer,edges_staight_inner,edges_staight_outer,edges_inside)
        return cell_edges_level_0,cell_edges_level_1,cell_edges_level_2

    cell_edges_level_0,cell_edges_level_1,cell_edges_level_2=calculateFreshPair(Names=Names,a1=a1,vertices_hex=vertices_hex,vertices_outer=vertices_outer,vertices_circle=vertices_circle,vertex_origin=vertex_origin)

    def partionByEdges(cellEdges,a2=a2,a1=a1,e1=e1,a=a,Names=Names):
        p1 = (0.0,a2,0.0)
        p2 = (-a1,a2,0.0) 
        midpoint = mid_point(p1,p2,edge=False)
        for pickedEdges in cellEdges:
            pickedCells = a.instances[Names.INSTANCE_NAME].cells
            a.PartitionCellByExtrudeEdge(line=e1.findAt(coordinates=midpoint), 
                                        cells=pickedCells, edges=pickedEdges, sense=const.REVERSE)

    partionByEdges(cellEdges=cell_edges_level_1)
    partionByEdges(cellEdges=cell_edges_level_2)
    partionByEdges(cellEdges=cell_edges_level_0)

def extrudePartionSketch2(Model,r_f,a1,a2,a3,Names):
    a = Model.rootAssembly
    r_i = r_f/2.0
    r_ii = r_i / 2.0
    def mid_point(p1,p2,a1=a1):
        if len(p1) == 2:
            p1 = (a1,p1[1],p1[0])
        if len(p2) == 2:
            p2 = (a1,p2[1],p2[0])
        return tuple([(x1+x2)/2.0 for x1,x2 in zip(p1,p2)])

    def circle_line_segment_intersection(pt1, pt2,circle_center, radius,  full_line=False, tangent_tol=1e-9):
        """ Find the points at which a circle intersects a line-segment.  This can happen at 0, 1, or 2 points.
        :param full_line: True to find intersections along full line - not just in the segment.  False will just return intersections within the segment.
        :param tangent_tol: Numerical tolerance at which we decide the intersections are close enough to consider it a tangent
        :return Sequence[Tuple[float, float]]: A list of length 0, 1, or 2, where each element is a point at which the circle intercepts a line segment.
        Note: We follow: http://mathworld.wolfram.com/Circle-LineIntersection.html
        """
        (p1x, p1y), (p2x, p2y), (cx, cy) = pt1, pt2, circle_center
        (x1, y1), (x2, y2) = (p1x - cx, p1y - cy), (p2x - cx, p2y - cy)
        dx, dy = (x2 - x1), (y2 - y1)
        dr = (dx ** 2 + dy ** 2)**.5
        big_d = x1 * y2 - x2 * y1
        discriminant = radius ** 2 * dr ** 2 - big_d ** 2
        if discriminant < 0:  # No intersection between circle and line
            return []
        else:  # There may be 0, 1, or 2 intersections with the segment
            intersections = [
                (cx + (big_d * dy + sign * (-1 if dy < 0 else 1) * dx * discriminant**.5) / dr ** 2,
                cy + (-big_d * dx + sign * abs(dy) * discriminant**.5) / dr ** 2)
                for sign in ((1, -1) if dy < 0 else (-1, 1))]  # This makes sure the order along the segment is correct
            if not full_line:  # If only considering the segment, filter out intersections that do not fall within the segment
                fraction_along_segment = [(xi - p1x) / dx if abs(dx) > abs(dy) else (yi - p1y) / dy for xi, yi in intersections]
                intersections = [pt for pt, frac in zip(intersections, fraction_along_segment) if 0 <= frac <= 1]
            if len(intersections) == 2 and abs(discriminant) <= tangent_tol:  # If line is tangent to circle, return just one point (as both intersections have same location)
                return [intersections[0]]
            else:
                return intersections
        
    def find_by_intersection(p1,p2,radius=r_f,z_c=a3,y_c=0.0):
        line_start = np.array([p1[1],p1[2]])
        line_end = np.array([p2[1],p2[2]])
        circle_center = np.array([y_c,z_c])
        l = circle_line_segment_intersection(line_start, line_end, circle_center=circle_center, radius=radius)
        return (p1[0],l[0][0],l[0][1])
    
    def calculateAllVertices(a1=a1,a2=a2,a3=a3,r_f=r_f):
        # the numbering is anticlockwise
        vertices_hex = []
        ang = math.cos(math.pi/4.0)
        vertices_hex.append((a1,0.0,a3-r_i*ang))
        vertices_hex.append((a1,r_i*ang,a3-r_i*ang))
        p1 = (a1,r_i*ang,a3-r_i*ang)
        p2 = (a1,r_i,a3)
        vertices_hex.append(mid_point(p1,p2))
        vertices_hex.append(p2)
        p1 = p2
        p2 = (a1,r_i*ang,a3+r_i*ang)
        vertices_hex.append(mid_point(p1,p2))
        vertices_hex.append(p2)
        vertices_hex.append((a1,0.0,a3+r_i*ang))
        p1 = (a1,-r_i*ang,a3+r_i*ang)
        vertices_hex.append(p1)
        p2 = (a1,-r_i,a3) 
        vertices_hex.append(mid_point(p1,p2))
        vertices_hex.append(p2)
        p1 = p2
        p2 = (a1,-r_i*ang,a3-r_i*ang)
        vertices_hex.append(mid_point(p1,p2))
        vertices_hex.append(p2)
        # vertices of the outer most geometry
        vertices_outer = []
        vertices_outer.append((a1,0.0,0.0))
        vertices_outer.append((a1,a2/2.0,0.0))
        vertices_outer.append((a1,a2,0.0))
        vertices_outer.append((a1,a2,a3))
        vertices_outer.append((a1,a2,2*a3))
        vertices_outer.append((a1,a2/2.0,2*a3))
        vertices_outer.append((a1,0.0,2*a3))
        vertices_outer.append((a1,-a2/2.0,2*a3))
        vertices_outer.append((a1,-a2,2*a3))
        vertices_outer.append((a1,-a2,a3))
        vertices_outer.append((a1,-a2,0.0))
        vertices_outer.append((a1,-a2/2.0,0.0))
        # vertices at the circle
        vertices_circle = []
        vertices_circle.append((a1,0.0,a3-r_f))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[1],p2=vertices_outer[1]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[2],p2=vertices_outer[2]))
        vertices_circle.append((a1,r_f,a3))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[4],p2=vertices_outer[4]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[5],p2=vertices_outer[5]))
        vertices_circle.append((a1,0.0,r_f + a3))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[7],p2=vertices_outer[7]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[8],p2=vertices_outer[8]))
        vertices_circle.append((a1,-r_f, a3))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[10],p2=vertices_outer[10]))
        vertices_circle.append(find_by_intersection(p1=vertices_hex[11],p2=vertices_outer[11]))
        # vertex origin
        vertex_origin = (a1,0.0,a3)
        vertices={}
        vertices['hex'] = vertices_hex
        vertices['outer'],vertices['circle'],vertices['origin'] = vertices_outer,vertices_circle,vertex_origin
        return vertices
    
    def arc_midpoint(start, end, center):
        if type(start) == tuple:
            start = np.array(list(start))
        if type(end) == tuple:
            end = np.array(list(end))
        if type(center) == tuple:
            center = np.array(list(center))
        # Calculate the vector from start to end point
        AB = end - start
        # Calculate the normal vector
        n = np.cross(AB, center - start)
        # Calculate the vector perpendicular to both AB and n
        v = np.cross(AB, n)
        # Normalize v to have length equal to the radius of the circular arc
        v = v / np.linalg.norm(v) * np.linalg.norm(center - start)
        # Calculate the midpoint
        midpoint = center + v
        return tuple(midpoint)
    
    vertices = calculateAllVertices(a1=a1,a2=a2,a3=a3,r_f=r_f)
    def angle2point(angle,radius,center=vertices['origin']):
        angle = angle * np.pi/180.0
        y = center[1] + radius * np.sin(angle)
        z = center[2] - radius * np.cos(angle)
        x = center[0]
        return (x,y,z)
    
    
    def facePartionCreation(coordinate,a1=a1,a2=a2,a=a,Names=Names):
        f1 = a.instances[Names.INSTANCE_NAME].faces
        e1 = a.instances[Names.INSTANCE_NAME].edges
        edgeIndices = f1.findAt(coordinate,).getEdges()
        pickedEdges = tuple([a.instances[Names.INSTANCE_NAME].edges[k] for k in edgeIndices])
        pickedCells = a.instances[Names.INSTANCE_NAME].cells
        p1 = (0.0,a2,0.0)
        p2 = (-a1,a2,0.0) 
        extrudeDirection = mid_point(p1,p2)
        a.PartitionCellByExtrudeEdge(line=e1.findAt(coordinates=extrudeDirection), 
                                    cells=pickedCells, edges=pickedEdges, sense=const.REVERSE)
    
    angles_innner = [45,135,225,315]
    for angle in angles_innner:
        facePartionCreation(angle2point(angle,r_ii))
    
    def find_edges_by_vertices(list1,list2,vertex_origin=vertices['origin'],a1=a1,Names=Names):
        list1_next = list1[1:]
        list1_next.append(list1[0])
        list2_next = list2[1:]
        list2_next.append(list2[0])
        list_midPoints = []
        def calAngle(p):
            return np.degrees(np.arctan2(p[1]-vertex_origin[1],-p[2]+vertex_origin[2]))
        for index,(p1,p2,p3,p4) in enumerate(zip(list1,list2,list1_next,list2_next)):
            midpoint1 = mid_point(p1,p2,a1=a1)
            midpoint2 = mid_point(p3,p4,a1=a1)
            angle1 = calAngle(midpoint1)
            angle2 = calAngle(midpoint2)
            if angle1<0: angle1 = 360 + angle1
            if angle2<0: angle2 = 360 + angle2
            if index == len(list1)-1: continue
            if index == len(list1)-1: angle2 = 360
            angle = (angle2-angle1)/2.0 + angle1 
            # print(index,' , angle1=',angle1,' ,angle2=',angle2,' ,angle=',angle)
            r1 = np.linalg.norm(np.array(list(map(lambda x,y:y-x,midpoint1,vertex_origin))))
            r2 = np.linalg.norm(np.array(list(map(lambda x,y:y-x,midpoint2,vertex_origin))))
            r = min(r1,r2)
            # print(index,' , r1=',r1,' ,r2=',r2,' ,r=',r)
            midpoint=angle2point(angle,r,center=vertex_origin)
            # a.ReferencePoint(point=midpoint)
            list_midPoints.append(midpoint)
            facePartionCreation(coordinate=midpoint)
            # print(index,' --- complete')
        return list_midPoints
    
    circle_mids = find_edges_by_vertices(list1=vertices['hex'],list2=vertices['circle'],vertex_origin=vertices['origin'],a1=a1,Names=Names)
    # i=0;k=f1.findAt(circle_mids[i],).index;highlight(f1[k:k+1])
    outer_faces = find_edges_by_vertices(list1=vertices['circle'],list2=vertices['outer'],vertex_origin=vertices['origin'],a1=a1,Names=Names)

def extrudePartionSketch3(vertices,Model,r_f,a1,a2,a3,Names):
    a = Model.rootAssembly
    r_i = r_f/2.0
    r_ii = r_i / 2.0
    def mid_point(p1,p2,a1=a1):
        if len(p1) == 2:
            p1 = (a1,p1[1],p1[0])
        if len(p2) == 2:
            p2 = (a1,p2[1],p2[0])
        return tuple([(x1+x2)/2.0 for x1,x2 in zip(p1,p2)])
    def angle2point(angle,radius,center=vertices['origin']):
        angle = angle * np.pi/180.0
        y = center[1] + radius * np.sin(angle)
        z = center[2] - radius * np.cos(angle)
        x = center[0]
        return (x,y,z)
    def facePartionCreation(coordinate,a1=a1,a2=a2,a=a,Names=Names):
        f1 = a.instances[Names.INSTANCE_NAME].faces
        e1 = a.instances[Names.INSTANCE_NAME].edges
        edgeIndices = f1.findAt(coordinate,).getEdges()
        pickedEdges = tuple([a.instances[Names.INSTANCE_NAME].edges[k] for k in edgeIndices])
        pickedCells = a.instances[Names.INSTANCE_NAME].cells
        p1 = (0.0,a2,0.0)
        p2 = (-a1,a2,0.0) 
        extrudeDirection = mid_point(p1,p2)
        a.PartitionCellByExtrudeEdge(line=e1.findAt(coordinates=extrudeDirection), 
                                    cells=pickedCells, edges=pickedEdges, sense=const.REVERSE)
    
    angles_innner = [45,135,225,315]
    for angle in angles_innner:
        facePartionCreation(angle2point(angle,r_ii))
    
    def find_edges_by_vertices(list1,list2,vertex_origin=vertices['origin'],a1=a1,Names=Names):
        list1_next = list1[1:]
        list1_next.append(list1[0])
        list2_next = list2[1:]
        list2_next.append(list2[0])
        list_midPoints = []
        def calAngle(p):
            return np.degrees(np.arctan2(p[1]-vertex_origin[1],-p[2]+vertex_origin[2]))
        for index,(p1,p2,p3,p4) in enumerate(zip(list1,list2,list1_next,list2_next)):
            midpoint1 = mid_point(p1,p2,a1=a1)
            midpoint2 = mid_point(p3,p4,a1=a1)
            angle1 = calAngle(midpoint1)
            angle2 = calAngle(midpoint2)
            if angle1<0: angle1 = 360 + angle1
            if angle2<0: angle2 = 360 + angle2
            if index == len(list1)-1: continue
            if index == len(list1)-1: angle2 = 360
            angle = (angle2-angle1)/2.0 + angle1 
            # print(index,' , angle1=',angle1,' ,angle2=',angle2,' ,angle=',angle)
            r1 = np.linalg.norm(np.array(list(map(lambda x,y:y-x,midpoint1,vertex_origin))))
            r2 = np.linalg.norm(np.array(list(map(lambda x,y:y-x,midpoint2,vertex_origin))))
            r = min(r1,r2)
            # print(index,' , r1=',r1,' ,r2=',r2,' ,r=',r)
            midpoint=angle2point(angle,r,center=vertex_origin)
            # a.ReferencePoint(point=midpoint)
            list_midPoints.append(midpoint)
            facePartionCreation(coordinate=midpoint)
            # print(index,' --- complete')
        return list_midPoints
    
    circle_mids = find_edges_by_vertices(list1=vertices['hex'],list2=vertices['circle'],vertex_origin=vertices['origin'],a1=a1,Names=Names)
    # i=0;k=f1.findAt(circle_mids[i],).index;highlight(f1[k:k+1])
    outer_faces = find_edges_by_vertices(list1=vertices['circle'],list2=vertices['outer'],vertex_origin=vertices['origin'],a1=a1,Names=Names)

def partitionGeometry(Model,r_f,a1,a2,a3,Names):
    depth = Names.GEOMETRY_DIVISION_DEPTH
    vertices = createSketchSingleChord( Model=Model,r_f=r_f,a1=a1,a2=a2,a3=a3,Names=Names,depth = depth)
    importSketchAndPartition(Model=Model,a1=a1,a2=a2,a3=a3,Names=Names)
    extrudePartionSketch3(   vertices=vertices,Model=Model,r_f=r_f,a1=a1,a2=a2,a3=a3,Names=Names)
