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

prefix = '20171229_QCF/1304eV_8e+16_ek'
names =['BC','QCF']

for name in names:
    datafile = os.path.join(prefix+ name,'processed_data/spectra/time-integrated_total')
    data = np.loadtxt(datafile,skiprows=1)
    plt.plot(data[:,0],data[:,1], label=name)


plt.legend()
plt.show()
