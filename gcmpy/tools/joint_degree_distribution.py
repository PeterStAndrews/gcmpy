
import networkx as nx

def joint_degree_distribution(G: nx.Graph):
    """
    Records the joint degree from a given network. Expects 
    the network nodes to have key 'joint_degree' which
    yields a tuple, this property is given when created using 
    the GCMAlgorithm class.
    :param nx.Graph:
    :return callable: degree distribution
    """
    num_vertices = G.order()
    
    PK_cache = {}
    for n in G.nodes():
        key = tuple(G.nodes[n]['joint_degree'])
        PK_cache[key] = PK_cache.get(key,0) + (1.0 / num_vertices)
  
    def p(key):
        return PK_cache.get(key,0)
    
    return p