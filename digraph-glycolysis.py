import os,sys,inspect
import KEGGutils as kg
from tabulate import tabulate
import networkx as nx
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

class Taxa:
    def __init__(self, list_V, list_E):
        self.Vertices = list_V
        self.Edges = list_E


def removeDuplicates(lst):
    return [t for t in (set(tuple(i) for i in lst))]

def generateEdgeList(speciesReaction):
    #storing substrate-product pairs (determining edges) in a list
    list_E = []

    for r_id, r_info in speciesReaction.items():
        #print("NodeID: ", r_id)
        temp_tuple = (r_info['substrates'][0][0], r_info['products'][0][0])
        #print(r_info['substrates'][0][0])
        #print(r_info['products'][0][0])

        #append tuple to edge list
        list_E.append(temp_tuple)

    #remove duplicate
    list_E = removeDuplicates(list_E)
    #print(list_E)

    return list_E

def generateVertexList(list_E):

    #create list of vertices
    list_V = [item for t in list_E for item in t]
    list_V = list(dict.fromkeys(list_V))
    #print(list_V)

    return list_V

def generateGraph(pathway):
    #draw graph
    #pathway.draw()

    listE = generateEdgeList(pathway.reactions)
    listV = generateVertexList(listE)
    speciesA = Taxa(listV, listE)

    return speciesA

def combinePathways(GA, GB):
    #combining list of edges/vertices from both Glycolysis and Citrate Cycle
    masterV = GA.Vertices + GB.Vertices
    #removing duplicates
    masterV = list(dict.fromkeys(masterV))
    masterE = removeDuplicates(GA.Edges + GB.Edges)
    GMaster = Taxa(masterV, masterE)

    return GMaster

def print_phenogram(title, dist_matrix, graph_names):
    clustered_matrix = hierarchy.linkage(dist_matrix, 'single')
    plt.figure()
    plt.title(title)
    pn = hierarchy.dendrogram(clustered_matrix, labels=graph_names)
    plt.show()

def bitCheck(DE, G):
    #DE is a (tuple) distinct edge from a list of tuples
    #G is a (class) distinct graph from a list of graphs
    if(DE in G.Edges):
        return 1
    else:
        return 0

def generateFPTable(G):
    #collate all distinct edges from the list of graphs
    FPTable = []
    DE = []

    for i in range(len(G)):
        for j in range(len(G[i].Edges)):
            DE.append(G[i].Edges[j])

    DE = removeDuplicates(DE)
    # print(DE)
    # print(type(DE))
    # print(len(DE))


    for x in range(len(DE)):
        bitDE = []
        row = []
        parity = 0
        row.append(DE[x])

        for y in range(len(G)):
            #check for presence of DE
            bitDE.append(bitCheck(DE[x], G[y]))

            #count parity
            parity = bitDE.count(1)

        row.append(bitDE)
        row.append(parity)
        FPTable.append(row)

    return FPTable


def sortTable(FPTable):
    #Sorting FP Table based on Parity
    FPTable.sort(key = lambda FPTable: FPTable[2], reverse=True)
    # print(FPTable)

    #Sorting FP Table based on decimal values
    for i in range(len(FPTable)):
        if (FPTable[i][2] != 5):
            # print(FPTable[i][1])
            #getting decimal values
            val = 0
            for j in range(len(FPTable[i][1])):
                val += FPTable[i][1][j] * (len(FPTable[i][1]) - j)
            FPTable[i][1].append(val)
            # print(FPTable[i][1])
        else:
            FPTable[i][1].append(15)
    FPTable.sort(key = lambda FPTable: FPTable[1][5], reverse=True)

    #remove unnecessary elements to display
    for i in range(len(FPTable)):
        FPTable[i].pop()
        FPTable[i][1].pop()
        FPTable[i].insert(0, i+1)

    return FPTable

def getInfoKegg(name, glycolysis_code, citrate_code):
    #getting info: Homo sapiens (human
    #selecting pathway: Glycolysis
    kg.get_infos(glycolysis_code)
    pathwayA = kg.KEGGpathway(pathway_id = glycolysis_code)
    GA = generateGraph(pathwayA)

    #selecting pathway: Citrate Cycle
    kg.get_infos(citrate_code)
    pathwayB = kg.KEGGpathway(pathway_id = citrate_code)
    GB = generateGraph(pathwayB)

    G = combinePathways(GA, GB)
    # print(f'Graph of {name}')
    # print(G.Vertices)
    # print(G.Edges)

    return G

def getJSISim(graphs, fptable):
    sim_mat = [[] for i in range(len(graphs))]
    for i in range(len(graphs)):
        for j in range(len(graphs)):
            curr_intersection = 0
            curr_union        = 0
            
            graph_1_code = ''
            graph_2_code = ''

            for row in fptable:
                bitcode = row[2]
                graph_1_code = graph_1_code + str(bitcode[i])
                graph_2_code = graph_2_code + str(bitcode[j])
            
            for k in range(len(graph_1_code)):
                if graph_1_code[k] == '1' or graph_2_code[k] == '1':
                    curr_union += 1

                if graph_1_code[k] == '1' and graph_2_code[k] == '1':
                    curr_intersection += 1

            sim_mat[i].append(curr_intersection / curr_union)

    print("Sim MAT")
    print(tabulate(sim_mat))

    return sim_mat

def getHammingSim(graphs, fptable):
    sim_mat = [[] for i in range(len(graphs))]

    for i in range(len(graphs)):
        for j in range(len(graphs)):
            curr_h_dist  = 0
            graph_1_code = ''
            graph_2_code = ''

            for row in fptable:
                bitcode = row[2]
                graph_1_code = graph_1_code + str(bitcode[i])
                graph_2_code = graph_2_code + str(bitcode[j])
            
            for k in range(len(graph_1_code)):
                if graph_1_code[k] != graph_2_code[k]:
                    curr_h_dist += 1

            sim_mat[i].append(curr_h_dist)

    return sim_mat

def main():
    #preliminaries
    currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    parentdir = os.path.dirname(currentdir)
    sys.path.insert(0,parentdir)

    print("Our default download directory is {}".format(kg.DOWNLOAD_DIR))
    ############

    #Set 1
    G1 = getInfoKegg("Homo sapiens", "hsa00010", "hsa00020")
    G2 = getInfoKegg("Schizosaccharomyces pombe", "spo00010", "spo00020")
    G3 = getInfoKegg("Ciona intestinalis (sea squirt)", "cin00010", "cin00020")
    G4 = getInfoKegg("Brassica oleracea", "boe00010", "boe00020")
    G5 = getInfoKegg("Enterobacter sichuanensis", "esh00010", "esh00020")


    ##### ANALYZING GRAPHS #####
    G_set1 = [G1,G2,G3,G4,G5]

    #Printing FP Table
    FPTable = generateFPTable(G_set1)
    # print(FPTable)

    #Sorting/Beautifying FP Table
    FPTable1_sorted = sortTable(FPTable)

    print(tabulate(FPTable1_sorted, headers=["SI No.", "DE", "BitCode"], tablefmt='grid'))


    #Set 2
    G6 = getInfoKegg("Vulcanisaeta distributa", "thg00010", "thg00020")
    G7 = getInfoKegg("Haemophilus influenzae", "hic00010", "hic00020")
    G8 = getInfoKegg("Ipomoea nil (Japanese morning glory)", "ini00010", "ini00020")
    G9 = getInfoKegg("Opisthorchis viverrini", "ovi00010", "ovi00020")
    G10 = getInfoKegg("Phytophthora sojae", "psoj00010", "psoj00020")


    ##### ANALYZING GRAPHS #####
    G_set2 = [G6,G7,G8,G9,G10]

    # for i in range(len(G_set2)):
    #     print(f"G{i+6} V: {len(G_set2[i].Vertices)} | G{i+6} E: {len(G_set2[i].Edges)}")


    #FP Table
    FPTable2 = generateFPTable(G_set2)
    FPTable2_sorted = sortTable(FPTable2)
    # print(tabulate(FPTable2_sorted, headers=["SI No.", "DE", "BitCode"], tablefmt='grid'))
    
    set1_names = ['hsa', 'spo', 'cin', 'boe', 'esh']
    set2_names = ['thg', 'hic', 'ini', 'ovi', 'psoj']

    set1_sim_mat_hd = getHammingSim(G_set1, FPTable1_sorted)
    set1_sim_mat_jsi = getJSISim(G_set1, FPTable1_sorted)
    set2_sim_mat_hd = getHammingSim(G_set2, FPTable2_sorted)
    set2_sim_mat_jsi = getJSISim(G_set2, FPTable2_sorted)

    print_phenogram("Set 1 - SIM using HD",set1_sim_mat_hd, set1_names)
    print_phenogram("Set 1 - SIM using JSI",set1_sim_mat_jsi, set1_names)
    
    print_phenogram("Set 2 - SIM using HD",set2_sim_mat_hd, set2_names)
    print_phenogram("Set 2 - SIM using JSI",set2_sim_mat_jsi, set2_names)


if __name__=="__main__":
    main()
