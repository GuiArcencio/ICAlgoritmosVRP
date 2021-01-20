import networkx as nx
import matplotlib.pyplot as plt
import sys

G = nx.DiGraph()

if len(sys.argv) > 1:
    points_file = sys.argv[1].strip()
    with open(points_file, 'r') as input_data_file:
        input_data = input_data_file.read()

    positions = {}
    demands = {}
    colors = []
    input_data = input_data.split('\n')
    nums = input_data[0].split(' ')
    N = int(nums[0])
    V = int(nums[1])
    C = float(nums[2])

    for i in range(1, N+1):
        nums = input_data[i].split(' ')
        positions[i-1] = (float(nums[1]), float(nums[2]))
        demands[i-1] = nums[0].strip()
        G.add_node(i-1)
        if i == 1: colors.append('#000000')
        else: colors.append('#59bdf0')

    emptyG = G.copy()

    with open(points_file + '.sol', 'r') as input_data_file:
        input_data = input_data_file.read()
    input_data = input_data.split('\n')
    
    # Edges
    for i in range(1, V + 1):
        clients = input_data[i].split(' ')
        if len(clients) > 2:
            for c in range(len(clients)-1):
                G.add_edge(int(clients[c]), int(clients[c+1]))

    plt.figure(1)
    nx.draw(emptyG, pos=positions, labels=demands, node_color = colors)
    plt.figure(2)
    nx.draw(G, pos=positions, labels=demands, node_color = colors)
    plt.show()

else:
    print('File path required.')