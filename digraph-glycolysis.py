import os,sys,inspect
import KEGGutils as kg
from tabulate import tabulate
import networkx as nx

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

    ##### [1] HOMO SAPIENS #######
    #getting info: Homo sapiens (human
    #selecting pathway: Glycolysis
    #kg.get_infos("hsa00010")
    pathway1A = kg.KEGGpathway(pathway_id = "hsa00010")
    G1A = generateGraph(pathway1A)

    #selecting pathway: Citrate Cycle
    #kg.get_infos("hsa00020")
    pathway1B = kg.KEGGpathway(pathway_id = "hsa00020")
    G1B = generateGraph(pathway1B)

    G1 = combinePathways(G1A, G1B)
    # print("G1 Master")
    # print(G1.Vertices)
    # print(G1.Edges)

    ##### end HOMO SAPIENS #######




    ##### [2] SCHIZOSACCHAROMYCES POMBE #######
    #getting info: Schizosaccharomyces pombe
    #selecting pathway: Glycolysis
    #kg.get_infos("spo00010")
    pathway2A = kg.KEGGpathway(pathway_id = "spo00010")
    G2A = generateGraph(pathway2A)

    #selecting pathway: Citrate Cycle
    #kg.get_infos("spo00020")
    pathway2B = kg.KEGGpathway(pathway_id = "spo00020")
    G2B = generateGraph(pathway2B)

    G2 = combinePathways(G2A, G2B)
    # print("G2 Master")
    # print(G2.Vertices)
    # print(G2.Edges)

    ##### end SCHIZOSACCHAROMYCES POMBE #######




    ##### [3] CIONA INTESTINALIS #######
    #getting info: Ciona intestinalis (sea squirt)
    #selecting pathway: Glycolysis
    #kg.get_infos("cin00010")
    pathway3A = kg.KEGGpathway(pathway_id = "cin00010")
    G3A = generateGraph(pathway3A)

    #selecting pathway: Citrate Cycle
    #kg.get_infos("cin00020")
    pathway3B = kg.KEGGpathway(pathway_id = "cin00020")
    G3B = generateGraph(pathway3B)

    G3 = combinePathways(G3A, G3B)
    # print("G3 Master")
    # print(G3.Vertices)
    # print(G3.Edges)

    ##### end SERINUS CANARIA #######




    ##### [4] BRASSICA OLERACEA #######
    #getting info: Brassica oleracea (wild cabbage)
    #selecting pathway: Glycolysis
    #kg.get_infos("boe00010")
    pathway4A = kg.KEGGpathway(pathway_id = "boe00010")
    G4A = generateGraph(pathway4A)
    # print("PRINTING HERE")
    # print(G4A.Vertices)
    # print(G4A.Edges)

    #selecting pathway: Citrate Cycle
    #kg.get_infos("boe00020")
    pathway4B = kg.KEGGpathway(pathway_id = "boe00020")
    G4B = generateGraph(pathway4B)

    G4 = combinePathways(G4A, G4B)
    # print("G4 Master")
    # print(G4.Vertices)
    # print(G4.Edges)

    ##### end AILUROPODA MELANOLEUCA #######




    ##### [5] ENTEROBACTER SICHUANENSIS #######
    #getting info: Enterobacter sichuanesis
    #selecting pathway: Glycolysis
    #kg.get_infos("esh00010")
    pathway5A = kg.KEGGpathway(pathway_id = "esh00010")
    G5A = generateGraph(pathway5A)

    #selecting pathway: Citrate Cycle
    #kg.get_infos("esh00020")
    pathway5B = kg.KEGGpathway(pathway_id = "esh00020")
    G5B = generateGraph(pathway5B)

    G5 = combinePathways(G5A, G5B)
    # print("G5 Master")
    # print(G5.Vertices)
    # print(G5.Edges)

    ##### end PHYSETER CATODON #######


    ##### ANALYZING GRAPHS #####
    G_set1 = [G1,G2,G3,G4,G5]

    # for i in range(len(G_set1)):
    #     print(f"G{i+1} V: {len(G_set1[i].Vertices)} | G{i+1} E: {len(G_set1[i].Edges)}")


    #Printing FP Table
    FPTable = generateFPTable(G_set1)
    print(FPTable)

    #Sorting/Beautifying FP Table
    FPTable1_sorted = sortTable(FPTable)

    print(tabulate(FPTable1_sorted, headers=["SI No.", "DE", "BitCode"], tablefmt='grid'))


    #Set 2
    ##### [6] Canis lupus familiaris #######
    G6 = getInfoKegg("dog", "cfa00010", "cfa00020")

    ##### [7] Octopus bimaculoides #######
    G7 = getInfoKegg("octopus", "obi00010", "obi00020")

    ##### [8] Delphinapterus leucas #######
    G8 = getInfoKegg("beluga whale", "dle00010", "dle00020")

    ##### [9] Panthera pardus #######
    G9 = getInfoKegg("leopard", "ppad00010", "ppad00020")

    ##### [10] Panthera tigris altaica #######
    G10 = getInfoKegg("amur tiger", "ptg00010", "ptg00020")

    ##### ANALYZING GRAPHS #####
    G_set2 = [G6,G7,G8,G9,G10]

    for i in range(len(G_set2)):
        print(f"G{i+6} V: {len(G_set2[i].Vertices)} | G{i+6} E: {len(G_set2[i].Edges)}")


    #FP Table
    FPTable2 = generateFPTable(G_set2)
    FPTable2_sorted = sortTable(FPTable2)
    print(tabulate(FPTable2_sorted, headers=["SI No.", "DE", "BitCode"], tablefmt='grid'))
    
    print(getHammingSim(G_set1, FPTable1_sorted))
    print(getJSISim(G_set1, FPTable1_sorted))
    print(getHammingSim(G_set2, FPTable2_sorted))
    print(getJSISim(G_set2, FPTable2_sorted))



if __name__=="__main__":
    main()
