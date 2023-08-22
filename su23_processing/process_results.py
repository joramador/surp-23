import matplotlib.pyplot as plt
import numpy as np

def comm_size(path):
    file = open(path, 'r')

    sizes = dict()

    for line in file:
        l = line.split("\t")
        cid = l[0]
        comm = l[1].strip("[]\n").split(",")
        c_size = len(comm)
        if c_size not in sizes:
            sizes[c_size] = 1
        else:
            count = sizes[c_size]
            sizes[c_size] = count + 1
    keys = list(sizes.keys())
    keys.sort()
    sorted_sizes = {i : sizes[i] for i in keys}

    return sorted_sizes



# key is node id
# value is list of communities that node belongs to (communities have this node as overlap)
def overlap(path):
    file = open(path, 'r')
    overlap = dict()

    for line in file:
        l = line.split("\t")
        cid = l[0]
        comm = l[1].strip("[]\n").split(",")
        for node in comm:
            node = int(node.strip())
            if node not in overlap:
                overlap[node] = [cid]
            else:
                overlap[node] = overlap[node] + [cid]
    return overlap

# key is the number of communities a node belongs to
# value is count of how many nodes have that same number of communities
def num_comm_node_belongs(overlap):
    num_belong = dict()
    for node in overlap:
        num_comm = len(overlap[node])
        if num_comm not in num_belong:
            num_belong[num_comm] = 1
        else:
            num_belong[num_comm] = num_belong[num_comm] + 1

    keys = list(num_belong.keys())
    keys.sort()
    sorted_num_belong = {i : num_belong[i] for i in keys}
    return sorted_num_belong

def overlap_figure(filename):
    overlap_0 = overlap("tiles/results/"+filename+"/t_0/strong-communities-0-t0.0")
    overlap_03 = overlap("tiles/results/"+filename+"/t_0.3/strong-communities-0-t0.3")
    overlap_05 = overlap("tiles/results/"+filename+"/t_0.5/strong-communities-0-t0.5")
    overlap_07 = overlap("tiles/results/"+filename+"/t_0.7/strong-communities-0-t0.7")
        
    plt.scatter(*zip(*num_comm_node_belongs(overlap_0).items()), label = "t = 0", alpha=0.4)
    plt.plot(*zip(*num_comm_node_belongs(overlap_0).items()),alpha=0.2, linestyle="-")

    plt.scatter(*zip(*num_comm_node_belongs(overlap_03).items()), label = "t = 0.3", alpha=0.4)
    plt.plot(*zip(*num_comm_node_belongs(overlap_03).items()), alpha=0.2, linestyle="--")

    plt.scatter(*zip(*num_comm_node_belongs(overlap_05).items()), label = "t = 0.5", alpha=0.4)
    plt.plot(*zip(*num_comm_node_belongs(overlap_05).items()),alpha=0.2,  linestyle="-.")

    plt.scatter(*zip(*num_comm_node_belongs(overlap_07).items()), label = "t = 0.7", alpha=0.4)
    plt.plot(*zip(*num_comm_node_belongs(overlap_07).items()),alpha=0.2,  linestyle=":")

    plt.xlabel("# of communties a node belongs to", fontsize="25")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("count", fontsize="25")
    plt.xticks(fontsize="25"), plt.yticks(fontsize="25")
    plt.legend(prop = {'size' : 30})
    plt.show()

def comm_size_figure(filename):
    size_0 = comm_size("tiles/results/"+filename+"/t_0/strong-communities-0-t0.0")
    #size_01 = comm_size("tiles/results/"+filename+"/t_0.1/strong-communities-0-t0.1")
    size_03 = comm_size("tiles/results/"+filename+"/t_0.3/strong-communities-0-t0.3")
    size_05 = comm_size("tiles/results/"+filename+"/t_0.5/strong-communities-0-t0.5")
    size_07 = comm_size("tiles/results/"+filename+"/t_0.7/strong-communities-0-t0.7")
    #size_09 = comm_size("tiles/results/"+filename+"/t_0.9/strong-communities-0-t0.9")
    #size_1 = comm_size("tiles/results/"+filename+"/t_1/strong-communities-0-t1.0")
        
    plt.scatter(*zip(*size_0.items()), label = "t = 0", alpha=0.4)
    plt.plot(*zip(*size_0.items()),alpha=0.2, linestyle="-")

    plt.scatter(*zip(*size_03.items()), label = "t = 0.3", alpha=0.4)
    plt.plot(*zip(*size_03.items()), alpha=0.2, linestyle="--")

    plt.scatter(*zip(*size_05.items()), label = "t = 0.5", alpha=0.4)
    plt.plot(*zip(*size_05.items()),alpha=0.2,  linestyle="-.")

    plt.scatter(*zip(*size_07.items()), label = "t = 0.7", alpha=0.4)
    plt.plot(*zip(*size_07.items()),alpha=0.2,  linestyle=":")

    plt.xlabel("community size", fontsize="25")
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel("count", fontsize="25")
    plt.xticks(fontsize="25"), plt.yticks(fontsize="25")
    plt.legend(prop = {'size' : 30})
    plt.show()

#comm_size_figure("lesmis")
#comm_size_figure("moreno_train") # nxt is asoiaf
#comm_size_figure("asoiaf")
#comm_size_figure("astro-ph")
comm_size_figure("cond-mat-2005")

#overlap_figure("lesmis")
#overlap_figure("moreno_train")
#overlap_figure("asoiaf")
#overlap_figure("astro-ph")
#overlap_figure("cond-mat-2005")