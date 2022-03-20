
import networkx as nx

class JointDegreeDistributionFromNetwork:
    @staticmethod   
    def get_joint_degree_distribution(G: nx.Graph) -> dict:
        """
        Records the joint degree from a given network. Expects 
        the network nodes to have key 'joint_degree' which
        yields a tuple, this property is given when created using 
        the GCMAlgorithm class.
        :param nx.Graph:
        :return dict: degree distribution
        """
        num_vertices = G.order()
    
        PK = {}
        for n in G.nodes():
            key = tuple(G.nodes[n]['joint_degree'])
            PK[key] = PK.get(key,0) + (1.0 / num_vertices)
  
        return PK