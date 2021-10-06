import bpy

# clear scene
for obj in bpy.context.scene.objects:
     obj.select_set(obj.type == "MESH")
bpy.ops.object.delete()
     
 
# make web mesh
vertices = [(0, 0, 0), (20, 20, 0), (-20, 20, 0)]
edges = [(0, 1), (1, 2), (2, 0)] 
#faces = [(0, 1, 2)]
faces = []

vertex_indices = {} # create dictionary to hold indices of the vertices in the vertices array
vertex_indices["anchor_threads"] = [0, 1, 2] # anchor vertices are the first three vertices in the vertices array

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
    
    frame_threads_vertices = [] # hold indices of frame threads
    
    # frame thread around vertex 1
    vertices.append(interpolate_edge_vertex(vertex_1, vertex_2, weight))
    vertices.append(interpolate_edge_vertex(vertex_1, vertex_3, weight))
    edges.append((num_vertices, num_vertices + 1))
    
    # frame thread around vertex 2
    vertices.append(interpolate_edge_vertex(vertex_2, vertex_1, weight))
    vertices.append(interpolate_edge_vertex(vertex_2, vertex_3, weight))
    edges.append((num_vertices + 2, num_vertices + 3))
    
    # frame thread around vertex 3
    vertices.append(interpolate_edge_vertex(vertex_3, vertex_1, weight))
    vertices.append(interpolate_edge_vertex(vertex_3, vertex_2, weight))
    edges.append((num_vertices + 4, num_vertices + 5))
    for i in range(0, 6):
        frame_threads_vertices.append(num_vertices + i)
        
    vertex_indices["frame_threads"] = frame_threads_vertices # add frame thread vertices to dictionary
    
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


# create spiderweb
frame_threads(vertices[0], vertices[1], vertices[2]) # add frame threads
anchor_vertices = vertex_indices["anchor_threads"]
triangle_center(vertices[anchor_vertices[0]], vertices[anchor_vertices[1]], vertices[anchor_vertices[2]])




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