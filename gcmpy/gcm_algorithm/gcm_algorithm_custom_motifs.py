import random
from itertools import chain, repeat, starmap

from gcmpy.gcm_algorithm.gcm_algorithm import GCMAlgorithm
from gcmpy.network.edge_list import LightWeightEdgeList
from gcmpy.names.gcm_algorithm_names import GCMAlgorithmNames


class GCMAlgorithmCustomMotifs(GCMAlgorithm):
    def __init__(self, params: dict):
        self._motif_indices: list = []
        try:
            self._motif_indices = params[GCMAlgorithmNames.MOTIF_INDICES]
            super().__init__(params)
        except Exception as e:
            raise (f"Error in {self.__class__.__name__}: {e}")

    def partition(self, lst: list, n: int) -> list:
        """
        Yield successive n-sized partitions from list.
        Silent if lst % n != 0.
        """
        return [lst[i : i + n] for i in range(0, len(lst), n)]

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

        # split the stub lists into partitions of equal
        # size to the number of that topology required to
        # construct the motif
        partitions: list[list] = []
        for i, k_list in enumerate(stubs):
            partitions.append(self.partition(k_list, self._motif_sizes[i]))

        # self._motif_vertices is a list of lists that contain
        # the indices of the joint degree slots required to construct
        # the motif.

        gen = self.infinite_sequence()

        # for each motif type
        for j, motif_indexes in enumerate(self._motif_indices):

            # for each motif of that type
            kk: int = motif_indexes[0]
            num_motifs = (0.0 + len(stubs[kk])) / self._motif_sizes[kk]
            for k in range(int(num_motifs)):

                vertices: list = []
                # for each orbit in the motif
                for index in motif_indexes:
                    vertices.append(partitions[index].pop())

                vertices = [item for sublist in vertices for item in sublist]

                # build the motif edges from the vertices
                es: list = self._build_functions[j](vertices)

                # get the motif id
                id = next(gen)
                EdgeList.motif_id.extend([id] * len(es))

                if len(es) == 2:
                    # if 2-clique tuple annoyingly unpacks ... so re-pack it
                    EdgeList.edge_list.extend([es])
                    EdgeList.topologies.extend([self._edge_names[j]()])

                else:
                    EdgeList.edge_list.extend(es)
                    EdgeList.topologies.extend(self._edge_names[j]())

        return EdgeList
