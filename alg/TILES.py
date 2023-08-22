"""
    Created on 11/feb/2015
    @author: Giulio Rossetti
"""
import networkx as nx
import gzip
import datetime
import time
from future.utils import iteritems
import os

import sys
if sys.version_info > (2, 7):
    from io import StringIO
    from queue import PriorityQueue
else:
    from cStringIO import StringIO
    from Queue import PriorityQueue


__author__ = "Giulio Rossetti"
__contact__ = "giulio.rossetti@gmail.com"
__website__ = "about.giuliorossetti.net"
__license__ = "BSD"


class TILES(object):
    """
        TILES
        Algorithm for evolutionary community discovery
    """

    # added threshold param - ja
    def __init__(self, filename=None, g=None, ttl=float('inf'), obs=7, path="", start=None, end=None, threshold=0.5):
        """
            Constructor
            :param g: networkx graph
            :param ttl: edge time to live (days)
            :param obs: observation window (days)
            :param path: Path where generate the results and find the edge file
            :param start: starting date
            :param end: ending date
            :param threshold: threshold for nodes to become core
        """
        self.path = path
        self.ttl = ttl
        self.cid = 0
        self.actual_slice = 0
        self.threshold = threshold # ADDED - JA

        # initialize graph
        if g is None:
            self.g = nx.Graph()
        else:
            self.g = g
        self.splits = None
        self.spl = StringIO()
        self.base = os.getcwd()
        # CHANGED to now include threshold value in name of file - JA
        self.status = open("%s/%s/extraction_status_t%.1f.txt" % (self.base, path, self.threshold), "w")
        self.removed = 0
        self.added = 0
        self.filename = filename
        self.start = start
        self.end = end
        self.obs = obs
        self.communities = {}

    def execute(self):
        """
            Execute TILES algorithm
        """
        self.status.write(u"Started! (%s) \n\n" % str(time.asctime(time.localtime(time.time()))))
        self.status.flush()
        qr = PriorityQueue()
        with open(self.filename, 'r') as f:
            first_line = f.readline()

        actual_time = datetime.datetime.fromtimestamp(float(first_line.split("\t")[2]))
        last_break = actual_time
        f.close()

        count = 0

        #################################################
        #                   Main Cycle                  #
        #################################################

        # while streaming source active
        f = open(self.filename)
        for l in f:
            l = l.split("\t")
            self.added += 1

            # get new interaction
            e = {}
            u = int(l[0])
            v = int(l[1])
            we = float(l[3]) # ADDED - JA
            dt = datetime.datetime.fromtimestamp(float(l[2]))

            ## uncomment below to include weight from line in source
            # where edge weight is a tuple, with the first being the number of
            # occurences for that edge and the second being the actual edge weight
            e['weight'] = (1, we) # ADDED - JA
            #e['weight'] = 1 # original
            e["u"] = l[0]
            e["v"] = l[1]
            # month = dt.month

            #############################################
            #               Observations                #
            #############################################

            gap = dt - last_break
            dif = gap.days

            if dif >= self.obs:
                last_break = dt
                self.added -= 1

                print("New slice. Starting Day: %s" % dt)

                self.status.write(u"Saving Slice %s: Starting %s ending %s - (%s)\n" %
                                  (self.actual_slice, actual_time, dt,
                                   str(time.asctime(time.localtime(time.time())))))

                self.status.write(u"Edge Added: %d\tEdge removed: %d\n" % (self.added, self.removed))
                self.added = 1
                self.removed = 0

                actual_time = dt
                self.status.flush()

                self.splits = gzip.open("%s/%s/splitting-%d.gz" % (self.base, self.path, self.actual_slice), "wt", 3)
                self.splits.write(self.spl.getvalue())
                self.splits.flush()
                self.splits.close()
                self.spl = StringIO()

                self.print_communities()
                self.status.write(
                    u"\nStarted Slice %s (%s)\n" % (self.actual_slice, str(datetime.datetime.now().time())))

            if u == v:
                continue

            # Check if edge removal is required
            if self.ttl != float('inf'):
                qr.put((dt, (int(e['u']), int(e['v']), int(e['weight'][0])))) # CHANGED - JA
                self.remove(dt, qr)

            # If node does not exist, then add the missing node
            if not self.g.has_node(u):
                self.g.add_node(u)
                self.g.nodes[u]['c_coms'] = {}  # central # chagned to nodes - ja

            if not self.g.has_node(v):
                self.g.add_node(v)
                self.g.nodes[v]['c_coms'] = {} # changed to nodes - ja
                
            # If edge does not exist, then add the edge
            # orignal dealt with occurances 
            # weight is not tuple of form - (occur, weight)
            if self.g.has_edge(u, v):
                w = self.g.adj[u][v]["weight"][0]
                self.g.adj[u][v]["weight"] = (w + e['weight'][0], e['weight'][1]) # changed - JA
                continue
            else:
                self.g.add_edge(u, v)
                self.g.adj[u][v]["weight"] = e['weight']

            u_n = list(self.g.neighbors(u))
            v_n = list(self.g.neighbors(v))

            #############################################
            #               Evolution                   #
            #############################################

            # new community of peripheral nodes (new nodes)
            if len(u_n) > 1 and len(v_n) > 1:
                common_neighbors = set(u_n) & set(v_n)
                self.common_neighbors_analysis(u, v, common_neighbors)

            count += 1

        #  Last writing
        self.status.write(u"Slice %s: Starting %s ending %s - (%s)\n" %
                          (self.actual_slice, actual_time, actual_time,
                           str(time.asctime(time.localtime(time.time())))))
        self.status.write(u"Edge Added: %d\tEdge removed: %d\n" % (self.added, self.removed))
        self.added = 0
        self.removed = 0

        self.print_communities()
        self.status.write(u"Finished! (%s)" % str(time.asctime(time.localtime(time.time()))))
        self.status.flush()
        self.status.close()

    @property
    def new_community_id(self):
        """
            Return a new community identifier
            :return: new community id
        """
        self.cid += 1
        self.communities[self.cid] = {}
        return self.cid

    def remove(self, actual_time, qr):
        """
            Edge removal procedure
            :param actual_time: timestamp of the last inserted edge
            :param qr: Priority Queue containing the edges to be removed ordered by their timestamps
        """

        # dict of communities to change
        coms_to_change = {}
        at = actual_time

        # main cycle on the removal queue
        if not qr.empty():

            # t is of form (timestamp, edge)
            t = qr.get()

            # get timestamp of when added to qr
            timestamp = t[0]

            # retrive (u, v, t)
            e = (t[1][0],  t[1][1], t[1][2])

            delta = at - timestamp
            displacement = delta.days

            if displacement < self.ttl:
                qr.put((timestamp, t[1]))

            # the time to live for the edge is expired
            else:
                while self.ttl <= displacement:

                    self.removed += 1
                    u = int(e[0])
                    v = int(e[1])
                    if self.g.has_edge(u, v):

                        # adjusted - JA
                        occur = self.g.adj,[u][v]["weight"][0] # added
                        w = self.g.adj[u][v]["weight"][1]

                        # decreasing link weight if greater than one
                        # (multiple occurrence of the edge: remove only the oldest)
                        if occur > 1:
                            # adjusted - JA
                            self.g.adj[u][v]["weight"][0] = occur - 1
                            e = (u, v, (occur-1, w))
                            qr.put((at, e))

                        else:
                            # u and v shared communities
                            if len(list(self.g.neighbors(u))) > 1 and len(list(self.g.neighbors(v))) > 1:
                                coms = set(self.g.nodes[u]['c_coms'].keys()) & set(self.g.nodes[v]['c_coms'].keys())  # changed to nodes - ja

                                for c in coms:
                                    if c not in coms_to_change:
                                        cn = set(self.g.neighbors(u)) & set(self.g.neighbors(v))
                                        coms_to_change[c] = [u, v]
                                        coms_to_change[c].extend(list(cn))
                                    else:
                                        cn = set(self.g.neighbors(u)) & set(self.g.neighbors(v))
                                        coms_to_change[c].extend(list(cn))
                                        coms_to_change[c].extend([u, v])
                                        ctc = set(coms_to_change[c])
                                        coms_to_change[c] = list(ctc)
                            else:
                                # if there are less than n (2) neighbors remove from community CHANGE HERE
                                if len(list(self.g.neighbors(u))) < 2:
                                    coms_u = [x for x in self.g.nodes[u]['c_coms'].keys()]  # changed to nodes - ja
                                    for cid in coms_u:
                                        self.remove_from_community(u, cid)

                                if len(list(self.g.neighbors(v))) < 2:
                                    coms_v = [x for x in self.g.nodes[v]['c_coms'].keys()]  # changed to nodes - ja
                                    for cid in coms_v:
                                        self.remove_from_community(v, cid)

                            self.g.remove_edge(u, v)

                    if not qr.empty():

                        t = qr.get()

                        timestamp = t[0]
                        delta = at - timestamp
                        displacement = delta.days

                        e = t[1]

        # update of shared communities
        self.update_shared_coms(coms_to_change)

    def update_shared_coms(self, coms_to_change):
        # update of shared communities
        print("update shard coms")
        for c in coms_to_change:
            if c not in self.communities:
                continue

            c_nodes = self.communities[c].keys()
            
            # more than 3 nodes in community CHANGE HERE
            if len(c_nodes) > 3:

                # subgraph: The graph structure cannot be changed but node/edge attributes can and are shared with the original graph.
                sub_c = self.g.subgraph(c_nodes)
                c_components = nx.number_connected_components(sub_c)

                # unbroken community -- everything is still connected as one
                if c_components == 1:
                    to_mod = sub_c.subgraph(coms_to_change[c])
                    self.modify_after_removal(to_mod, c)

                # broken community: bigger one maintains the id, the others obtain a new one

                # many components exist -- biggest one is the actual community
                # remaining components are new nodes / components being introduced
                else:
                    new_ids = []

                    first = True

                    # sets of nodes, one for each component of G.
                    components = nx.connected_components(sub_c)
                    # for each set
                    for com in components:
                        if first:
                            # less than 3 -> broken community CHANGE HERE
                            print("broken community")
                            if len(com) < 3:
                                self.destroy_community(c)
                            else:
                                to_mod = list(set(com) & set(coms_to_change[c]))
                                sub_c = self.g.subgraph(to_mod)
                                self.modify_after_removal(sub_c, c)
                            first = False

                        else:
                            # more than 3 -> okay CHANGE HERE
                            if len(com) > 3:
                                # update the memberships: remove the old ones and add the new one
                                to_mod = list(set(com) & set(coms_to_change[c]))
                                sub_c = self.g.subgraph(to_mod)
                                print("more than 3")
                                central = self.centrality_test(sub_c).keys()
                                # 3? not sure about this CHANGE HERE
                                if len(central) >= 3:
                                    actual_id = self.new_community_id
                                    new_ids.append(actual_id)
                                    
                                    # n is a key
                                    #self.g.edges(central) #returns edges involved in core
                                    for n in central:
                                        # retrieve edges of node n
                                        #self.g.edges(n)
                                        self.add_to_community(n, actual_id)

                    # splits
                    if len(new_ids) > 0 and self.actual_slice > 0:
                        self.spl.write(u"%s\t%s\n" % (c, str(new_ids)))
            else:
                self.destroy_community(c)

    def modify_after_removal(self, sub_c, c):
        """
            Maintain the clustering coefficient invariant after the edge removal phase
            :param sub_c: sub-community to evaluate
            :param c: community id
        """
        central = self.centrality_test(sub_c).keys()

        # in case of previous splits, update for the actual nodes
        remove_node = set(self.communities[c].keys()) - set(sub_c.nodes())

        for rm in remove_node:
            self.remove_from_community(rm, c)

        # removing community 3 CHANGE HERE
        if len(central) < 3:
            self.destroy_community(c)
        else:
            not_central = set(sub_c.nodes()) - set(central)
            for n in not_central:
                self.remove_from_community(n, c)

    def common_neighbors_analysis(self, u, v, common_neighbors):
        """
            General case in which both the nodes are already present in the net.
            :param u: a node
            :param v: a node
            :param common_neighbors: common neighbors of the two nodes
        """

        # no shared neighbors
        if len(common_neighbors) < 1:
            return

        else:

            shared_coms = set(self.g.nodes[v]['c_coms'].keys()) & set(self.g.nodes[u]['c_coms'].keys()) # changed to nodes - ja
            only_u = set(self.g.nodes[u]['c_coms'].keys()) - set(self.g.nodes[v]['c_coms'].keys())  # changed to nodes - ja
            only_v = set(self.g.nodes[v]['c_coms'].keys()) - set(self.g.nodes[u]['c_coms'].keys())  # changed to nodes - ja

            # community propagation: a community is propagated iff at least two of [u, v, z] are central
            propagated = False
            print("core candidate involving nodes u:" + str(u) + " v: " +str(v) )
            for z in common_neighbors:
                # ADDED WEIGHT CALCULATIONS - JA
                uv_w = self.g[u][v]["weight"][1]
                uz_w = self.g[u][z]["weight"][1]
                vz_w = self.g[v][z]["weight"][1]
                hrmnc = (3*uv_w*uz_w*vz_w)/((uv_w*uz_w)+(uv_w*vz_w)+(uz_w*vz_w))
                normalized_hrmnc = hrmnc / max(uv_w, uz_w, vz_w)
                #print("normalized hrmnc: " + str(normalized_hrmnc))
                if (normalized_hrmnc >= self.threshold):
                    for c in self.g.nodes[z]['c_coms'].keys():  # changed to nodes - ja
                    
                        if c in only_v:
                            # u is not part of community c
                            self.add_to_community(u, c)
                            propagated = True

                        if c in only_u:
                            # v is not part of community c
                            self.add_to_community(v, c)
                            propagated = True

                    for c in shared_coms:
                        # z is not part of community c
                        if c not in self.g.nodes[z]['c_coms']:  # changed to nodes - ja
                            self.add_to_community(z, c)
                            propagated = True

            else:
                if not propagated:
                    # new community
                    uv_w = self.g[u][v]["weight"][1]
                    #print("weight of edge" +str(u)+str(v)+ " is " + str(uv_w))
                    for z in common_neighbors:
                        uz_w = self.g[u][z]["weight"][1]
                        #print("weight of edge" +str(u)+str(z)+ " is " + str(uz_w))
                        vz_w = self.g[v][z]["weight"][1]
                        #print("weight of edge" +str(v)+str(z)+ " is " + str(vz_w))
                        hrmnc = (3*uv_w*uz_w*vz_w)/((uv_w*uz_w)+(uv_w*vz_w)+(uz_w*vz_w))
                        normalized_hrmnc = hrmnc / max(uv_w, uz_w, vz_w)
                        # if threshold is met, make it core
                        #print("normalized hrmnc: " + str(normalized_hrmnc))
                        if (normalized_hrmnc >= self.threshold):
                            actual_cid = self.new_community_id
                            self.add_to_community(u, actual_cid)
                            self.add_to_community(v, actual_cid)
                            self.add_to_community(z, actual_cid)

    def print_communities(self):
        """
            Print the actual communities
        """
        # CHANGED to now include threshold value in name of file - JA
        out_file_coms = gzip.open("%s/%s/strong-communities-%d-t%.1f.gz" % (self.base, self.path, self.actual_slice, self.threshold), "wt", 3)
        com_string = StringIO()

        nodes_to_coms = {}
        merge = {}
        coms_to_remove = []
        drop_c = []

        self.status.write(u"Writing Communities (%s)\n" % str(time.asctime(time.localtime(time.time()))))
        self.status.flush()
        for idc, comk in iteritems(self.communities):

            com = comk.keys()

            if self.communities[idc] is not None:
                # more than 2 okay CHANGE HERE
                if len(com) > 2:
                    key = tuple(sorted(com))

                    # Collision check and merge index build (maintaining the lowest id)
                    if key not in nodes_to_coms:
                        nodes_to_coms[key] = idc
                    else:
                        old_id = nodes_to_coms[key]
                        drop = idc
                        if idc < old_id:
                            drop = old_id
                            nodes_to_coms[key] = idc

                        # merged to remove
                        coms_to_remove.append(drop)
                        if not nodes_to_coms[key] in merge:
                            merge[nodes_to_coms[key]] = [idc]
                        else:
                            merge[nodes_to_coms[key]].append(idc)
                else:
                    drop_c.append(idc)
            else:
                drop_c.append(idc)

        write_count = 0
        for k, idk in iteritems(nodes_to_coms):
            write_count += 1
            if write_count % 50000 == 0:
                out_file_coms.write(com_string.getvalue())
                out_file_coms.flush()
                com_string = StringIO()
                write_count = 0
            com_string.write(u"%d\t%s\n" % (idk, str(list(k))))

        for dc in drop_c:
            self.destroy_community(dc)

        out_file_coms.write(com_string.getvalue())
        out_file_coms.flush()
        out_file_coms.close()

        # write the graph
        self.status.write(u"Writing actual graph status (%s)\n" % str(time.asctime(time.localtime(time.time()))))
        self.status.flush()
        # CHANGED to now include threshold value in name of file - JA
        out_file_graph = gzip.open("%s/%s/graph-%d-t%.1f.gz" % (self.base, self.path, self.actual_slice, self.threshold), "wt", 3)
        g_string = StringIO()
        for e in self.g.edges():
            g_string.write(u"%d\t%s\t%d\t%d\n" % (e[0], e[1], self.g.adj[e[0]][e[1]]['weight'][0], self.g.adj[e[0]][e[1]]['weight'][1])) # changed to include occurance and weight - ja

        out_file_graph.write(g_string.getvalue())
        out_file_graph.flush()
        out_file_graph.close()

        # Write merge status
        self.status.write(u"Writing merging file (%s)\n" % str(time.asctime(time.localtime(time.time()))))
        self.status.flush()
        # CHANGED to now include threshold value in name of file - JA
        out_file_merge = gzip.open("%s/%s/merging-%d-t%.1f.gz" % (self.base, self.path, self.actual_slice, self.threshold), "wt", 3)
        m_string = StringIO()
        for comid, c_val in iteritems(merge):
            # maintain minimum community after merge
            c_val.append(comid)
            k = min(c_val)
            c_val.remove(k)
            m_string.write(u"%d\t%s\n" % (k, str(c_val)))
        out_file_merge.write(m_string.getvalue())
        out_file_merge.flush()
        out_file_merge.close()

        # Community Cleaning
        m = 0
        for c in coms_to_remove:
            self.destroy_community(c)
            m += 1

        self.status.write(u"Merged communities: %d (%s)\n" % (m, str(time.asctime(time.localtime(time.time())))))

        self.actual_slice += 1
        self.status.write(u"Total Communities %d (%s)\n" % (len(self.communities.keys()),
                                                           str(time.asctime(time.localtime(time.time())))))
        self.status.flush()

    def destroy_community(self, cid):
        nodes = [x for x in self.communities[cid].keys()]
        for n in nodes:
            self.remove_from_community(n, cid)
        self.communities.pop(cid, None)

    # ADD THE CONDITION OF MEETING THRESHOLD WHERE THIS FUNCTION APPEARS
    def add_to_community(self, node, cid):

        # assign cid to the node in graph
        self.g.nodes[node]['c_coms'][cid] = None  # changed to nodes - ja

        # cid in dict of communities (i.e. community exsists)
        if cid in self.communities:
            # add node as a key to the community
            self.communities[cid][node] = None
        else:
            # make cid as a key and its value be a 
            # new dict for the node in its community
            self.communities[cid] = {node: None}

    def remove_from_community(self, node, cid):
        if cid in self.g.nodes[node]['c_coms']:  # changed to nodes - ja
            self.g.nodes[node]['c_coms'].pop(cid, None)  # changed to nodes - ja
            if cid in self.communities and node in self.communities[cid]:
                self.communities[cid].pop(node, None)

    def centrality_test(self, subgraph):
        print("centrality test")
        central = {}

        for u in subgraph.nodes():

            # if node u not in core yet
            if u not in central:
                cflag = False
                neighbors_u = set(self.g.neighbors(u))

                # if node u has more than one neighbor,
                # meaning that it is connected to two other nodes
                if len(neighbors_u) > 1:
                    for v in neighbors_u:
                        # v is a neighbor of u
                        # if the two nodes are not the same
                        if u > v:

                            # if core, then we are done
                            if cflag:
                                break
                            else:
                                neighbors_v = set(self.g.neighbors(v))
                                # cn = neighbors_v & neighbors_v # ORIGINAL BUT SHOULDNT THIS BE U?
                                cn = neighbors_u & neighbors_v
                                if len(cn) > 0:

                                    # ORIGINAL
                                    # central[u] = None
                                    # central[v] = None
                                    uv_w = self.g[u][v]["weight"][1]

                                    for n in cn:
                                        # ADDED THESE LINES - JA
                                        # get edge weight u,n
                                        un_w = self.g[u][n]["weight"][1]
                                        # get edge weight v,n
                                        vn_w = self.g[v][n]["weight"][1]

                                        # calculating value
                                        hrmnc = (3*uv_w*un_w*vn_w)/((uv_w*un_w)+(uv_w*vn_w)+(un_w*vn_w))
                                        normalized_hrmnc = hrmnc / max(uv_w, un_w, vn_w)
                                        # if threshold is met, make it core
                                        print("normalized hrmnc: " + str(normalized_hrmnc))
                                        if (normalized_hrmnc > self.threshold):
                                            central[u] = None
                                            central[v] = None
                                            central[n] = None
                                        # ORIGINAL
                                        # central[n] = None
                                    # core is established
                                    cflag = True
        return central


# hacky way to run file: change filename here to run this modified TILES algorithm
an = TILES(filename="tiles/datasets/astro-ph/astro-ph.tsv", threshold=0)
an.execute()
print("number of unique cid's "+str(an.cid))
