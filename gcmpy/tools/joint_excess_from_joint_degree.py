
from gcmpy.tools.average_joint_degree_from_jdd import AverageJointDegreeFromJDD

class JointExcessDDfromJDD:
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
                    q[tuple(_joint_degree)] = ((
                        (_joint_degree[index] + 1) * jdd[joint_degree] + 0.0) / averages[index]
                    )
            qks.append(q)
        return qks




