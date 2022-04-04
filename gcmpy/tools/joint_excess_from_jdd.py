from gcmpy.tools.average_joint_degree_from_jdd import AverageJointDegreeFromJDD


class JointExcessfromJDD:
    @staticmethod
    def get_joint_excess_distributions(jdd: dict) -> list[dict]:

        qks = []

        averages = AverageJointDegreeFromJDD.get_average_joint_degrees(jdd)

        joint_degrees = list(jdd.keys())
        num_topologies: int = len(joint_degrees[0])
        for index in range(num_topologies):
            q = {}
            for joint_degree in joint_degrees:
                _joint_degree = list(joint_degree)
                if _joint_degree[index] > 0:
                    _joint_degree[index] -= 1
                    q[tuple(_joint_degree)] = (
                        (_joint_degree[index] + 1) * jdd[joint_degree] + 0.0
                    ) / averages[index]
            qks.append(q)
        return qks

    @staticmethod
    def convert_list_qks_to_dict(qks_list: list[dict], keys: list[str]) -> dict[dict]:
        """
        Static method to convert a list of dicts to dict of dicts
        :param qks_list: list of excess degree distributions
        :param keys: list of topologies
        :returns dict: converted dict object
        """
        qks_dict = {}
        for key, qk in zip(keys, qks_list):
            qks_dict[key] = qk
        return qks_dict

    @staticmethod
    def convert_dict_qks_to_list(qks_dict: dict[dict], keys: list[str]) -> list[dict]:
        """
        Static method to convert a dict of dicts to list of dicts
        :param qks_dict: dict keyed by topology string
        :param keys list: list of topology keys (ordered)
        :returns list[dict]:
        """
        qks_list: list = []
        for key in keys:
            qks_list.append(qks_dict[key])
        return qks_list
