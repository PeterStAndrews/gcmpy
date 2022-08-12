class JointDegreeFromExcess:
    @staticmethod
    def invert_single(qk: dict, i: int) -> dict:
        """
        Invert a excess joint degree distribution to a joint
        degree distribution.
        :param qk: dict excess joint degree distribution
        :param i: index of current topology in joint degree tuple
        :returns dict: joint degree distribution
        """
        P = {}
        bottom = sum(
            [(qk[joint_excess] / (joint_excess[i] + 1)) for joint_excess in qk]
        )
        for joint_excess in qk:
            top = qk[joint_excess] / (joint_excess[i] + 1)
            joint_degree = list(joint_excess)
            joint_degree[i] += 1
            P[tuple(joint_degree)] = top / bottom
        return P

    @staticmethod
    def observations_from_dict(qks: dict, keys: list) -> dict:
        P_observations = {}
        for i, key in enumerate(keys):
            P_observations[key] = JointDegreeFromExcess.invert_single(qks[key], i)
        return P_observations

    @staticmethod
    def get_joint_degree_distribution(qks: dict, keys: list) -> dict:
        p_obs: dict = JointDegreeFromExcess.observations_from_dict(qks, keys)
        # get set of common keys across all observations
        p_obs_list = [p_obs[topology] for topology in p_obs]
        common_keys = list(set.intersection(*map(set, p_obs_list)))

        if not common_keys:
            raise "Error: JointDDFromExcess - no common keys found across observations"

        # choose a common key and a reference topology to scale to
        common_key = common_keys[0]
        choesn_topology = "2-clique"

        # scale all observations to the chosen value
        base_value = p_obs[choesn_topology][common_key]
        for topology in p_obs:
            if choesn_topology == topology:
                continue
            scale_factor = base_value / p_obs[topology][common_key]
            for key in p_obs[topology]:
                p_obs[topology][key] *= scale_factor

        # merge keys
        P = {}
        for topology in p_obs:
            P.update(p_obs[topology])

        # renormalise
        total = sum(P.values())
        for k in P:
            P[k] /= total

        return P
