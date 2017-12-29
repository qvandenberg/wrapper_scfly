#__________________________________________________________________________
#
#  Created by Q.Y. van den Berg on 07/11/2017.
#
#  Process SCFLY output with "new_scfly_analyse.py" in conjunction with "new_lib.analyse.py" Library
#  Use this script to plot results
#
#  python3 plot_scfly.py "output directory"
#__________________________________________________________________________

import sys, os, getopt
import matplotlib.pyplot as plt
import numpy as np

weights_Q = np.loadtxt('fscanweights',skiprows=1)
weights_S = np.loadtxt('fscanweights_Sam',skiprows=1)

x_q = weights_Q[:,1]
y_q =weights_Q[:,2]
x_s = weights_S[:,1]
y_s =weights_S[:,2]

plt.plot(x_q,y_q, label="Quincy")
plt.plot(x_s,y_s,label="Sam")
plt.legend()
plt.show()
