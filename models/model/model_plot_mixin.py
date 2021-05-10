import matplotlib.pyplot as plt
import numpy as np
from models.model.domain import Domains
from models.model.variable import Variables


class ModelPlotMixin:

    domains = Domains()
    variables = Variables()

    def plot_result(self, t, y):

        n_domains = len(list(self.domains.keys()))

        if n_domains == 1:

            x_key = list(self.domains.keys())[0]

            y0 = y[0]

            x = self.domains[x_key]

            for i, ti in enumerate(t):
                if i > 0:
                    yi = y[i]
                    keys = list(self.variables.keys())
                    fig, axs = plt.subplots(len(keys))
                    for j, var in enumerate(keys):
                        if len(keys) == 1:
                            axsj = axs
                        else:
                            axsj = axs[j]

                        yj = self.variables[var].parse(yi)
                        y0j = self.variables[var].parse(y0)
                        axsj.plot(x, y0j, "k-")
                        axsj.plot(x, yj, "bo-")
                        axsj.set_ylabel(var)
                        axsj.set_xlabel(x_key)
                        axsj.legend(["t={}".format(t[0]), "t={}".format(t[i])])

        elif n_domains == 2:

            x_keys = list(self.domains.keys())

            x1 = self.domains[x_keys[0]]
            x2 = self.domains[x_keys[1]]
            X1, X2 = np.meshgrid(x1, x2, indexing="ij")

            for i, ti in enumerate(t):
                if i > 0:
                    Yi = np.reshape(y[i], (len(x1), len(x2)))
                    keys = list(self.variables.keys())
                    fig, axs = plt.subplots(subplot_kw={"projection": "3d"})
                    for j, var in enumerate(keys):
                        if len(keys) == 1:
                            axsj = axs
                        else:
                            axsj = axs[j]
                        axsj.plot_surface(X1, X2, Yi)
                        axsj.set_xlabel(x_keys[0])
                        axsj.set_ylabel(x_keys[1])
                        axsj.set_zlabel(var)
                    plt.title('For t={}s'.format(ti))

        else:

            print("Exceeds expected number of domains")

        plt.show()
