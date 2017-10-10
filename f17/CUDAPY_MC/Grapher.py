from matplotlib import pyplot
from AOMC import *
from Ihelp import *

instance = Run()
paths = instance.ree()
path_plot_range = min(Number_of_Paths, 50)
for i in range(path_plot_range):
    pyplot.plot(paths[i])
pyplot.show()

