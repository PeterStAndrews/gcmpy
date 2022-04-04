import random
from iteration_utilities import grouper
from itertools import chain, repeat, starmap

from gcmpy.gcm_algorithm.gcm_algorithm import GCMAlgorithm
from gcmpy.network.edge_list import LightWeightEdgeList


class GCMAlgorithmFast(GCMAlgorithm):
    def random_clustered_graph(self, jds: list) -> LightWeightEdgeList:

        stubs = [
            list(chain.from_iterable(starmap(repeat, r)))
            for r in map(enumerate, zip(*jds))
        ]

        # shuffle each stub list
        for k_list in stubs:
            random.shuffle(k_list)

        # create list for edges and add joint degree sequence
        EdgeList = LightWeightEdgeList()
        EdgeList.joint_degrees = jds

        gen = self.infinite_sequence()

        # for each topology list ...
        for k, k_list in enumerate(stubs):

            # iterate the degree list
            for vertices in grouper(k_list, self._motif_sizes[k]):
                # add the edges to the network using the builder callback
                es = self._build_functions[k](list(vertices))
                EdgeList.edge_list.extend(es)

                # add the edge names to a list
                EdgeList.topologies.extend([self._edge_names[k]] * len(es))

                # record the motif id
                id = next(gen)
                EdgeList.motif_id.extend([id] * len(es))

        # return the graph model as a LightWeightEdgeList
        return EdgeList
