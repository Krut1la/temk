class Interpolator(object):
    def __init__(self, x_nodes, y_nodes):
        """
        :param x_nodes: x values in nodes
        :param y_nodes: y values in nodes
        """
        self._x_nodes = x_nodes
        self._y_nodes = y_nodes

    @classmethod
    def get_method_name(cls):
        pass

    def interpolate(self, x):
        """
        Evaluates polynomial in x point
        :param x: point to evaluate
        :return: y value
        """
        pass

    def generate_y_range(self, x_range):
        """
        """
        y = []
        for x in x_range:
            y.append(self.interpolate(x))

        return y


class InterpolatorNewton(Interpolator):
    def __init__(self, x_nodes, y_nodes):
        super(InterpolatorNewton, self).__init__(x_nodes, y_nodes)
        self.method_name = "Newton"
        self.__end_diffs = self.__get_poly_newton_end_diffs()

    @classmethod
    def get_method_name(cls):
        return "Newton"

    def __get_poly_newton_end_diffs(self):
        """
        Cals end diffs of Newton polynomial using the general formula (1.6)
        :return: array of coefficients
        """
        m = len(self._x_nodes)

        end_diffs = []

        for n in range(1, m + 1):
            S = 0
            for j in range(1, n + 1):
                P = 1
                for i in range(1, n + 1):
                    if i == j:
                        continue

                    P = P * (self._x_nodes[j - 1] - self._x_nodes[i - 1])

                S = S + self._y_nodes[j - 1] / P

            end_diffs.append(S)

        return end_diffs

    def interpolate(self, x):
        n = len(self._x_nodes) - 1
        p = self.__end_diffs[n]
        for k in range(1, n + 1):
            p = self.__end_diffs[n - k] + (x - self._x_nodes[n - k]) * p

        return p