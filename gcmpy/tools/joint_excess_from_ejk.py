

class JointExcessDDFromEjk:
    """
    Extracts the joint excess degree distributions from a 
    list of joint excess joint degree distributions. This is
    done by summing the rows of the ejk matrices.
    """
    @staticmethod
    def get_excess_joint_distributions(
        ejks: list[dict], excess_degree_keys: list[list]
    ) -> list[dict]:

        if len(ejks) != len(excess_degree_keys):
            raise('Error in JointExcessDDFromEjk: incorrect ejk or keys provided')

        qks: list = []
        for ejk, keys in zip(ejks,excess_degree_keys):
            q = {}
            for left_key in keys:
                for right_key in keys:
                    q[left_key] = q.get(left_key, 0.0) + ejk[left_key+right_key]
            qks.append(q)
        return qks