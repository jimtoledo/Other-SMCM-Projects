#Jim Toledo Graph Theory Project 2 - Minimum Roman Domination for Randomly Generated Graphs
import math
import networkx as nx
import matplotlib.pyplot as plt
import random

def connected(G, size, root): #DFS for all connected nodes in the graph
    visited[root] = True #set current node to true
    for i in range(size):
        if not visited[i] and G[root][i] == 1: #go through all unvisited nodes connected to the current node
            connected(G,size,i)

def isConnected(G, size, root): #returns whether or not the graph G is connected
    global visited
    visited = [False for i in range(size)] #initially set all nodes as unvisited
    connected(G, size, root) #call DFS
    return visited == [True for i in range(size)] #graph is connected if all nodes have been visited by the DFS

def genGraph(n,r): #generate random connected graph of order 'n' and sparsity around 'r'
    matrix = [[0 for i in range(n)] for j in range(n)] #initial empty (no edge) graph of order n
    for i in range(n):
        for j in range(i+1,n): #only iterate through one side of the diagonal of the adjacency matrix since it is symmetric
            if random.uniform(0,1) < r: #chance that any two vertices are connected is r
                matrix[i][j]=1
                matrix[j][i]=1
    while not isConnected(matrix, size, 0): #add random edges until graph is connected
        i = random.randint(0,n-1)
        j = random.randint(0,n-1)
        if i==j: #no self-loops
            continue
        matrix[i][j]=1
        matrix[j][i]=1
    
    return matrix

def isRomanDomination(G, zeros, ones, twos): #check if graph function f : V(G) -> {0, 1, 2} is a Roman dominating function
    for zero_vertex in zeros:
        connected_to_two = False
        for two_vertex in twos:
            if G[zero_vertex][two_vertex] == 1: #pretty much just check that all zero vertices are connected to a two vertex
                connected_to_two = True
                break
        if not connected_to_two:
            return False
    return True

def partition(G, V, num_zeros, num_ones, num_twos, zero_list, one_list, two_list): #generate and test all partitions of the vertex set V given #s needed in each set
    if len(V) > 0:
        vertex = V[0]
        #check all ways in which the vertex is in bin 0, bin 1, and bin 2
        if num_zeros > 0:
            zero_list.append(vertex) #color vertex 0 (red)
            copyV = V.copy()
            copyV.remove(vertex) #remove the vertex for recursive call so that it cannot be recolored
            part = partition(G, copyV, num_zeros - 1, num_ones, num_twos, zero_list, one_list, two_list)
            if part[0]:
                return (True, part[1], part[2], part[3]) #exit out of loop early when RD is found
            zero_list.remove(vertex) #if all combos where the vertex was colored '0' fails, remove coloring and try another
        if num_ones > 0:
            one_list.append(vertex)
            copyV = V.copy()
            copyV.remove(vertex)
            part = partition(G, copyV, num_zeros, num_ones - 1, num_twos, zero_list, one_list, two_list)
            if part[0]:
                return (True, part[1], part[2], part[3]) #exit out of loop early when RD is found
            one_list.remove(vertex)
            
        if num_twos > 0:
            two_list.append(vertex)
            copyV = V.copy()
            copyV.remove(vertex)
            part = partition(G, copyV, num_zeros, num_ones, num_twos - 1, zero_list, one_list, two_list)
            if part[0]:
                return (True, part[1], part[2], part[3]) #exit out of loop early when RD is found
            two_list.remove(vertex)
        return (False, [], [], []) #return false if no partition works
    else: #if no more vertices to color
        if isRomanDomination(G, zero_list, one_list, two_list):
            return (True, zero_list.copy(), one_list.copy(), two_list.copy()) #return the colorings and also the fact that the colorings produce an RD
        else:
            return (False, [], [], [])

def checkPartitionRD(G, num_zeros, num_ones, num_twos): #generate and test all partitions of the vertex set V given #s needed in each set
    V = list(range(len(G))) #vertex set of G
    return partition(G, V, num_zeros, num_ones, num_twos, [], [], [])

def getCombos(n): #all possible sizes of coloring sets given a graph of order n, in order from smallest to largest sum
    combos = []
    for num_twos in range(math.ceil(n/2)): #minimum RD can never contain any more than n/2 2's
        for num_ones in range(n - num_twos + 1):
            num_zeros = n - num_twos - num_ones
            if num_zeros == 0 or num_twos > 0: #impossible to have 0's and no 2's in the coloring
                combos.append((num_ones + 2 * num_twos, num_zeros, num_ones, num_twos)) #tuple where first element is the sum of vertices after coloring
    combos.sort() #test from lowest number to highest
    return combos

def findMinRD(G): #finds the minimum Roman Domination of graph G
    combos = getCombos(len(G)) #get all possible sizes of coloring sets
    for c in combos:
        pot_RD = checkPartitionRD(G, c[1], c[2], c[3])
        if pot_RD[0]: #if minimum RD is found
            return (pot_RD[1], pot_RD[2], pot_RD[3], c[0]) #return tuple of list of 0s, 1s, and 2s; the last element of the tuple is the RD number
    return ([], [], [], 0) #if no RD can be found

def convertList(l, dict): #used to label nodes for visual
    newList = []
    for e in l:
        newList.append(dict[e])
    return newList

def makeNetworkGraph(G): #make the NetworkX Graph from adjacency matrix representation
    n = len(G)
    G_vis = nx.Graph()
    G_vis.add_nodes_from(list(range(n)))
    for i in range(n):
        for j in range(n):
            if G[i][j] == 1: G_vis.add_edge(i, j)
    return G_vis
    
size = int(input("Order: "))
if size < 4: #no trivial RD problems >:(
    quit()
visited = [False for i in range(size)]

G = genGraph(size,0.3) #generate graph (adjacency matrix) from the user inputted size

G_vis = makeNetworkGraph(G)

labels = {}
for i in range(size):
    labels[i] = chr(i + 97) #label nodes a,b,c,...

print("Graph generated, now finding min. Roman Domination")

#Original Example [Britain, Gaul, Iberia, Rome, Constantinople, North Africa, Egypt, Asia Minor]
#G = [[0, 1, 1, 0, 0, 0, 0, 0],
#     [1, 0, 1, 1, 0, 0, 0, 0],
#     [1, 1, 0, 1, 0, 1, 0, 0],
#     [0, 1, 1, 0, 1, 1, 1, 0],
#     [0, 0, 0, 1, 0, 0, 1, 1],
#     [0, 0, 1, 1, 0, 0, 1, 0],
#     [0, 0, 0, 1, 1, 1, 0, 1],
#     [0, 0, 0, 0, 1, 0, 1, 0]]
#G_vis = makeNetworkGraph(G)
#labels[0]='B'
#labels[1]='G'
#labels[2]='I'
#labels[3]='R'
#labels[4]='C'
#labels[5]='N'
#labels[6]='E'
#labels[7]='A'

result = findMinRD(G)
zeros = convertList(result[0], labels)
ones = convertList(result[1], labels)
twos = convertList(result[2], labels)
print("0s (red):", *zeros, sep = " ")
print("1s (green):", *ones, sep = " ")
print("2s (blue):", *twos, sep = " ")
print("Roman Domination Number: ", str(result[3]))

pos = nx.spring_layout(G_vis)
nx.draw_networkx_nodes(G_vis, pos, result[0], 500, 'r')
nx.draw_networkx_nodes(G_vis, pos, result[1], 500, 'g')
nx.draw_networkx_nodes(G_vis, pos, result[2], 500, 'b')

nx.draw_networkx_labels(G_vis, pos, labels, font_size = 16)

nx.draw_networkx_edges(G_vis, pos, width = 1.0)

plt.axis('off')
plt.show()
