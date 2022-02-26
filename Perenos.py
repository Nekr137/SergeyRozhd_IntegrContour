from numpy import mean

class Perenos:
    def __init__(self, xPos) -> None:
        self._xPos = xPos
        self._depth_values = []

        # cached
        self._sorted_by_depths = []
        self._has_sorted_by_depths = False

    def append(self, depth, value):
        self._depth_values.append((depth, value))
        self._has_sorted_by_depths = False

    def get_sorted_by_depths(self):
        if not self._has_sorted_by_depths:
            self._sorted_by_depths = sorted(self._depth_values, key=lambda it: it[0])
            self._has_sorted_by_depths = True
        return self._sorted_by_depths

    def get_depths(self):
        return [v[0] for v in self._depth_values]
    
    def get_values(self):
        return [v[1] for v in self._depth_values]
    
    def get_min_depth(self):
        return self.get_sorted_by_depths()[0]

    def get_max_depth(self):
        return self.get_sorted_by_depths()[-1]

    def get_mean_values(self):
        return mean([v[1] for v in self._depth_values])

    def get_depth_step(self):
        dv = self.get_sorted_by_depths()
        depth_deriv = [dv[i+1][0] - dv[i][0] for i in range(len(dv) - 1)]
        depth_step = mean(depth_deriv)
        return depth_step
