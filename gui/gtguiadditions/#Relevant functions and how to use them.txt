"""
#Relevant functions and how to use them
from sage.all import Graph

#defining a graph
G = Graph() #defines the variable as a graph
H = Graph({1:[2,3], 2:[3]}) #creates graph with vertices 1,2,3 and edges 1-2, 1-3, 2-3

#adding vertices
G.add_vertex(1) #adds vertex named 1 to the graph G
G.add_vertex('a') #adds vertex named a to the graph G
G.add_vertices([1, 2, 3]) #can add multiple vertices at one time
G.add_vertices(['a', 'b', 'c'])

#remove vertex
G.delete_vertex(1) #deletes vertex one and all incident edges
G.delete_vertices(1, 2) #deletes multiple vertices and all incident edges

#adding edges
G.add_edges([(1,2)]) #add single edge from 1 to 2
G.add_edges([(1,2), (2,3)]) #adds edges from 1 to 2 and 2 to 3

#misc
G.show() #creates png of graph in seperate window

G.density() #shows density of the graph

G.degree() #displays the degrees as an array [1,2,1] etc

G.is_eulerian() #returns true or false
G.eulerian_circuit(labels=False) #shows path for ciruit [(1,2), (2,3)] etc

G.is_hamiltonian #returns true or false
G.hamiltonian_cycle() #need to look more in to. Does something weird

G.girth() #returns girth of the graph ex: 3

G.is_planar() #returns true or false

G.is_connected() #returns true or false
"""

