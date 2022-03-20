
import networkx as nx

class JointExcessJointDegreeDistribution:
    
    _topology_names: list[str]
    _G: nx.Graph = None
        
    _num_edges: list[int] = []                      # number of edges of i-th topology
    _ejks: list[dict] = []                          # list of ejk dicts/matrices
    _degree_keys: list = []                         # list of joint degree tuples
    _excess_degree_keys: list[list] = []            # list of excess degree tuple lists

    @property
    def excess_degree_keys(self) -> dict:
        return self._excess_degree_keys
    @excess_degree_keys.setter
    def excess_degree_keys(self, value: dict) -> None:
        self._excess_degree_keys = value

    def __init__(self, params: dict):
        try:
            self._G = params['network']
            self._topology_names: list[str] = params['edge_names']
            self._num_edges = [0] * len(params['edge_names'])
        
            if 'jdd' in params:
                self._degree_keys = params['jdd'].keys()
            else:
                self._degree_keys = set(
                    [tuple(self._G.nodes[n]['joint_degree']) for n in self._G.nodes()]
                )
        except Exception as e:
            raise (f"Error in JointExcessJointDegreeDistribution: {e}")

        self.resolve_excess_degree_keys()
    
    def resolve_excess_degree_keys(self) -> None:
        """
        Obtains the excess degree keys from the joint degree keys 
        for each topology.
        """
        self._excess_degree_keys = []
        for i, _ in enumerate(self._topology_names):
            excess_keys_top_i = []
            for jd in self._degree_keys:
                if jd[i] > 0:
                    jd = list(jd)
                    jd[i] -= 1
                    excess_keys_top_i.append(tuple(jd))
            self._excess_degree_keys.append(set(excess_keys_top_i))
            
    def count_edge_types(self) -> None:
        """
        Enumerates the number of edges of the i-th topology. 
        Assumes that the network edges have attribute `topology`.
        """
        for i, name in enumerate(self._topology_names):
            for e in self._G.edges():
                if self._G.edges[e]['topology'] == name:
                    self._num_edges[i] += 1
        
    def get_ejk(self, i: int, name: str) -> dict:
        """
        Get mixing patterns for the i-th topology. Routine assumes that 
        the network edges have attribute `topology` and that the vertices 
        have attribute `joint_degree`. 
        :param i int: index of current topology in joint degree
        :param name str: key for edge name
        :returns ejk dict: mixing patterns for this topology
        """
        ejk = {}
        for e in self._G.edges():
            if self._G.edges[e]['topology'] == name:
            
                u,v = e
            
                u_joint_degree = list(self._G.nodes[u]['joint_degree'])
                v_joint_degree = list(self._G.nodes[v]['joint_degree'])
            
                u_joint_degree[i] -= 1
                u_joint_excess_degree = tuple(u_joint_degree)
            
                v_joint_degree[i] -=1
                v_joint_excess_degree = tuple(v_joint_degree)
            
                key1 = u_joint_excess_degree + v_joint_excess_degree
                key2 = v_joint_excess_degree + u_joint_excess_degree
            
                # divide by 2*num_edges due to adding each one twice from each end
                ejk[key1] = ejk.get(key1,0) + (1.0 / (2*self._num_edges[i]))
                ejk[key2] = ejk.get(key2,0) + (1.0 / (2*self._num_edges[i]))
        
        return ejk
            
    def get_ejks(self) -> list[dict]:
        """
        Get mixing patterns for all edge topologies
        """
        self.count_edge_types()
        for i, name in enumerate(self._topology_names):
            self._ejks.append(self.get_ejk(i,name))
        return self._ejks
            