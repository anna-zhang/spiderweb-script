import bpy
from math import sqrt, pi, sin, cos
import numpy as np



# clear scene
#for obj in bpy.context.scene.objects:
#     obj.select_set(obj.type == "MESH")
#bpy.ops.object.delete()
     
 
# make web mesh
vertices = [(5, 8, 0), (9, 4, 0), (-2, 6, 0)] # need to error check that all three vertices aren't on a single line
edges = []
#edges = []
#faces = [(0, 1, 2)]
faces = []

vertex_indices = {} # create dictionary to hold indices of the vertices in the vertices array
vertex_indices["anchor_threads"] = [0, 1, 2] # anchor vertices are the first three vertices in the vertices array
edge_indices = {} # create dictionary to hold indices of the edges in the edges array


# returns the coordinates of the vertex along the edge connecting vertex_1 and vertex_2, according to a weight value (0.5 weight value gets halfway point)
def interpolate_edge_vertex(vertex_1, vertex_2, weight) -> ():
    vertex_x = vertex_1[0] + (vertex_2[0] - vertex_1[0]) * weight 
    vertex_y = vertex_1[1] + (vertex_2[1] - vertex_1[1]) * weight 
    vertex_z = vertex_1[2] + (vertex_2[2] - vertex_1[2]) * weight 
    return (vertex_x, vertex_y, vertex_z)


#print(interpolate_edge_vertex(vertices[0], vertices[1], 0.5))

# returns the coordinate of the center of a triangle formed by vertex_1, vertex_2, and vertex_3
def triangle_center(vertex_1, vertex_2, vertex_3):
    vertex_x = (vertex_1[0] + vertex_2[0] + vertex_3[0]) / 3.0
    vertex_y = (vertex_1[1] + vertex_2[1] + vertex_3[1]) / 3.0
    vertex_z = (vertex_1[2] + vertex_2[2] + vertex_3[2]) / 3.0
    
    vertices.append((vertex_x, vertex_y, vertex_z)) # add vertex for triangle center
    vertex_indices["center"] = len(vertices) - 1

#print(triangle_center(vertices[0], vertices[1], vertices[2]))

# takes in three vertices forming a triangle and gets six vertices along the edges of the triangle and sets those as the vertices forming the frame threads
def frame_threads(vertex_1, vertex_2, vertex_3):
    weight = 0.25 # go a quarter distance out from the first vertex along the edge toward the second vertex
    num_vertices = len(vertices) # get current number of vertices in web mesh to make sure we're joining the right vertices for edges
    num_edges = len(edges) # get current number of edges in web mesh
    
    frame_threads_vertices = [] # hold indices of frame threads vertices
    frame_threads_edges = [] # hold indices of frame threads edges
    
    # frame thread around vertex 1
    vertices.append(interpolate_edge_vertex(vertex_1, vertex_2, weight))
    vertices.append(interpolate_edge_vertex(vertex_1, vertex_3, weight))
    edges.append((num_vertices, num_vertices + 1)) # add frame thread connecting two frame vertices around vertex 1
    
    # frame thread around vertex 2
    vertices.append(interpolate_edge_vertex(vertex_2, vertex_1, weight))
    edges.append((num_vertices, num_vertices + 2)) # add frame thread along anchor vertex 1 and anchor vertex 2
    vertices.append(interpolate_edge_vertex(vertex_2, vertex_3, weight))
    edges.append((num_vertices + 2, num_vertices + 3)) # add frame thread connecting two frame vertices around vertex 2
    
    # frame thread around vertex 3
    vertices.append(interpolate_edge_vertex(vertex_3, vertex_1, weight))
    edges.append((num_vertices + 1, num_vertices + 4)) # add frame thread along anchor vertex 1 and anchor vertex 3
    vertices.append(interpolate_edge_vertex(vertex_3, vertex_2, weight))
    edges.append((num_vertices + 4, num_vertices + 5)) # add frame thread connecting two frame vertices around vertex 3
    edges.append((num_vertices + 3, num_vertices + 5)) # add frame thread along anchor vertex 2 and anchor vertex 3

    for i in range(0, 6):
        frame_threads_vertices.append(num_vertices + i) # save indices of added frame vertices
        frame_threads_edges.append(num_edges + i) # save indices of added frame edges
        
    vertex_indices["frame_threads"] = frame_threads_vertices # add frame thread vertices to dictionary
    edge_indices["frame_threads"] = frame_threads_edges # add frame thread edges to dictionary
    
    for i in range(0, 3): # add anchor edges (edges between anchor vertices and frame vertices)
        edges.append((vertex_indices["anchor_threads"][i], frame_threads_vertices[2 * i]))
        edges.append((vertex_indices["anchor_threads"][i], frame_threads_vertices[2 * i + 1]))
    edge_indices["anchor_threads"] = list(range(len(edges) - 1, len(edges) + 5)) # save indices of added anchor edges
    
    
#    # frame thread around vertex 1
#    vertices.append(interpolate_edge_vertex(vertex_1, vertex_2, weight))
#    vertices.append(interpolate_edge_vertex(vertex_1, vertex_3, weight))
#    edges.append(num_vertices, num_vertices + 1)
#    
#    # frame thread around vertex 2
#    vertices.append(interpolate_edge_vertex(vertex_2, vertex_1, weight))
#    vertices.append(interpolate_edge_vertex(vertex_2, vertex_3, weight))
#    # frame thread around vertex 3
#    vertices.append(interpolate_edge_vertex(vertex_3, vertex_1, weight))
#    vertices.append(interpolate_edge_vertex(vertex_3, vertex_2, weight))

def get_distance(v1, v2): # takes in two points and returns the distance between two points
    distance = sqrt((v1[0] - v2[0]) * (v1[0] - v2[0]) + (v1[1] - v2[1]) * (v1[1] - v2[1]) + (v1[2] - v2[2]) * (v1[2] - v2[2]))
    return distance

def triangle_area(vertex_1, vertex_2, vertex_3): # calculates the area of the triangle formed by vertex_1, vertex_2, and vertex_3
    vector_1 = np.asarray(tuple(map(lambda i, j: i - j, vertex_2, vertex_1))) # vertex_2 - vertex_1
    vector_2 = np.asarray(tuple(map(lambda i, j: i - j, vertex_3, vertex_1))) # vertex_2 - vertex_1
    magnitude = np.linalg.norm(np.cross(vector_1, vector_2)) # get magnitude of the cross product of vector_1 and vector_2
    area = 0.5 * magnitude
    return area

def outer_circle_starting_radius(): # compute the outer circle radius, largest circle centered at the triangle center such that it doesn't hit any of the frame or anchor edges
    # compute distance from triangle center to the closest frame vertex and use that as the starting outer radius
    center = vertices[vertex_indices["center"]]
    frame_vertices = vertex_indices["frame_threads"] # get all of the frame thread vertices
    shortest_distance = float('inf')
    for index in range(len(frame_vertices)):
        frame_vertex = vertices[vertex_indices["frame_threads"][index]] # get vertex location of the frame vertex
        distance = get_distance(center, frame_vertex) # compute distance from center to that frame vertex
        if distance < shortest_distance:
            shortest_distance = distance # update shortest distance if the distance between the center to that frame vertex is shorter than what's currently the shortest
    return shortest_distance


def circle_line_intersect(circle, v1, v2): # compute whether a line formed by vertex_1 and vertex_2 intersects a sphere with cross-section circle; return true if it hits the circle, false otherwise
    r = circle[3] # get radius, circle format (center_x, center_y, center_z, radius)
    print("checking radius r: " + str(r))
    a = ((v2[0] - v1[0]) * (v2[0] - v1[0])) + ((v2[1] - v1[1]) * (v2[1] - v1[1])) + ((v2[2] - v1[2]) * (v2[2] - v1[2])) 
    b = 2 * ((v2[0] - v1[0]) * (v1[0] - circle[0]) + (v2[1] - v1[1]) * (v1[1] - circle[1]) + (v2[2] - v1[2]) * (v1[2] - circle[2]))
    c = (circle[0] * circle[0]) + (circle[1] * circle[1]) + (circle[2] * circle[2]) + (v1[0] * v1[0]) + (v1[1] * v1[1]) + (v1[2] * v1[2]) - 2 * (circle[0] * v1[0] + circle[1] * v1[1] + circle[2] * v1[2]) - r * r
    if ((b * b) - (4 * a * c)) < 0:
        return False # line does not intersect circle
    else:
        return True # line intersects circle
    
    
def find_outer_circle_radius(center, radius): # computes the outer circle radius by taking the center of the triangle and a starting radius (computed from outer_circle_starting_radius)
    x, y, z = center # get the x, y, z of the triangle center vertex
    indices = edge_indices["frame_threads"] # get the indices of the edges of the frame threads in the edges array
    intersect = False # hold whether a circle centered at center with this radius intersects the frame threads
    for i in indices: # compute whether the circle intersects with the line (each frame thread edge)
        vertex_1_index, vertex_2_index = edges[i] # get the indices of the vertices forming the frame thread edge
        if (circle_line_intersect((x, y, z, radius), vertices[vertex_1_index], vertices[vertex_2_index])): # check intersection
            intersect = True # circle hits frame thread edge
            break
    if (intersect):
        # recursively call with a smaller radius
        return find_outer_circle_radius(center, radius / 1.2) # need to return for both cases when recursively calling function, or else will get back None
    else:
        return radius # found outer circle radius value
    

def create_circle(center, radius, num_vertices): # add vertices of circle with center "center" and radius "radius"
    theta = (2 * pi) / num_vertices # calculate angle of each slice of the circle
    circle_center = np.asarray(center)
    anchor_vertices = vertex_indices["anchor_threads"] # get indices of anchor vertices 
    vector_1 = np.asarray(vertices[anchor_vertices[1]]) - np.asarray(vertices[anchor_vertices[0]]) # vector 1 defining plane containing circle
    vector_2 = np.asarray(vertices[anchor_vertices[2]]) - np.asarray(vertices[anchor_vertices[0]]) # vector 2 defining plane containing circle
    v1v2_cross = np.cross(vector_1, vector_2) # get cross product of two vectors that form the triangle edges to get a vector normal to the plane with the circle
    v1v2_cross_normalized = v1v2_cross / np.linalg.norm(v1v2_cross) # normalize the vector that's normal to the plane containing the circle
    u = vector_1 / np.linalg.norm(vector_1) # get unit vector to serve as "x-axis" of the plane
    v = np.cross(u, v1v2_cross_normalized) # get vector to serve as "y-axis" of the plane
    v = v / np.linalg.norm(v) # normalize v to get the unit "y-axis" of the plane
    print("dot product: " + str(np.dot(u, v))) # make sure dot product of u and v = 0 so the axes are perpendicular
    num_vertices_total = len(vertices) # get current # of vertices
    circle_vertices = [] # save indices of the vertices lying on the circle
    for i in range(0, num_vertices):
        vertex = tuple(circle_center + radius * cos(theta * i) * u + radius * sin(theta * i) * v) 
        vertices.append(vertex)
        print("circle_vertex: " + str(vertex))
        circle_vertices.append(num_vertices_total + i) # save index in vertices list of the vertex being added
    vertex_indices["outer_circle"] = circle_vertices
    print(str(vertex_indices["outer_circle"]))
        
 
def connect_frame_and_outer_circle(): # create an edge between each frame vertex and its closest vertex on the outer circle
    frame_vertices = vertex_indices["frame_threads"] # get all of the indices of the frame thread vertices
    outer_circle_vertices = vertex_indices["outer_circle"] # get all of the indices of the outer circle vertices
    frame_outer_circle_edges = []
    total_edges = len(edges)
    
    for index in range(len(frame_vertices)): # iterate through every frame vertex
        frame_vertex_index = vertex_indices["frame_threads"][index] # get the index of the frame vertex in the vertices list
        frame_vertex = vertices[frame_vertex_index] # get vertex location of the frame vertex
        shortest_distance = float('inf')
        closest_circle_vertex_index = 0
        for i in range(len(outer_circle_vertices)): # iterate through every outer circle vertex
            outer_circle_vertex = vertices[vertex_indices["outer_circle"][i]] # get vertex location of the outer circle vertex 
            distance = get_distance(frame_vertex, outer_circle_vertex) # compute distance from frame vertex to outer circle vertex
            if distance < shortest_distance:
                shortest_distance = distance # update shortest distance if the distance between the frame vertex to that outer circle vertex is shorter than what's currently the shortest
                closest_circle_vertex_index = vertex_indices["outer_circle"][i] # update index to be that outer circle vertex's index
        # create edge between that frame vertex and its closest circle vertex
        edges.append([frame_vertex_index, closest_circle_vertex_index])
        frame_outer_circle_edges.append(total_edges + index) # save the index of the added edge in the edges list
    
    edge_indices["frame_outer_circle"] = frame_outer_circle_edges # save the list of indices of the edges between the frame vertices and their cloest outer circle vertex

# create spiderweb
frame_threads(vertices[0], vertices[1], vertices[2]) # add frame threads
anchor_vertices = vertex_indices["anchor_threads"] # get indices of anchor vertices 
triangle_center(vertices[anchor_vertices[0]], vertices[anchor_vertices[1]], vertices[anchor_vertices[2]]) # compute triangle center from the three anchor vertices
print("area: " + str(triangle_area(vertices[anchor_vertices[0]], vertices[anchor_vertices[1]], vertices[anchor_vertices[2]]))) # compute triangle area
print("outer_circle_starting_radius: " + str(outer_circle_starting_radius())) # compute outer circle starting radius
center_vertex = vertices[vertex_indices["center"]] # get circle center vertex
outer_radius = find_outer_circle_radius(center_vertex, outer_circle_starting_radius()) # compute outer circle radius
print("found outer_circle_radius: " + str(outer_radius))
create_circle(center_vertex, outer_radius, 8)
connect_frame_and_outer_circle()

#def threads_from_center():
#    num_vertices = len(vertices) # get current number of vertices in web mesh to make sure we're joining the right vertices for threads
#    vertices.append(triangle_center(vertices[0], vertices[1], vertices[2])) # add vertex for triangle center
#    center_vertex_index = num_vertices
#    starting_index = num_vertices - 6 # six vertices and index starts at 0
#    for i in range(0, 6):
#        edges.append((center_vertex_index, starting_index + i))
#        
#threads_from_center()




    
web_mesh = bpy.data.meshes.new('web_mesh')
web_mesh.from_pydata(vertices, edges, faces)
web_mesh.update()
# make web object from web mesh
web_object = bpy.data.objects.new('web_object', web_mesh)
# make web collection
web_collection = bpy.data.collections.new('web_collection')
bpy.context.scene.collection.children.link(web_collection)
# add web object to scene collection
web_collection.objects.link(web_object)