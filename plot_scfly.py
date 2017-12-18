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
import matplotlib.plot as plt
import numpy
from new_lib_analyse import *

## Check number of input arguments. Should be only 1.
if len(sys.argv) != 2:
    print( 'Number of arguments is ', len(sys.argv)-1, '; one required!.')

# Build folders for processed data
init = initFolderStructure(sys.argv[1])
input_parameters = loadInput(sys.argv[1])

## Time integrate emission spectra
spec = timeIntSpec(input_parameters) # (input parameters, i_start, i_end, t_start, t_end) only "input parameters" is mandatory
spec.time_integrate(input_parameters)

## Apply f-scan. Build in synchronisation/lock threads. Make sure spectra are written out before executing this part.
supergauss_parameters = [0.957467, 3.54691, 0.46181, 0.042533, 0.1929, 0.21648] #
spec.fscan(supergauss_parameters)

## Smoothening of spectra
# spec.broaden(np.linspace(1,20,20), 200, 'GAUSS') # (vector of folders, width [ev], lineshape)
# spec.broaden('TOTAL', 200, 'GAUSS') # (vector of folders, width [ev], lineshape)

## Extract: temperature-density conditions, population, rates
extract = extract(input_parameters)
extract.temperature_density(input_parameters)
extract.populations(input_parameters,[3,9],'gs')
extract.rates(input_parameters,[3,9],'gs','coll_ion')
# x=extract.superconfiguration(input_parameters.Z,[4,7],'single_ch')