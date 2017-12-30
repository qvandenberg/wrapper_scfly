#__________________________________________________________________________
#
#  Created by Q.Y. van den Berg on 07/11/2017.
#
#  Process SCFLY output to plot spectra, rates, population, temperature/density
#  python3 scfly_analyse.py "output directory"
#__________________________________________________________________________

import sys, os, getopt
import numpy
import re
from new_lib_analyse import *
from subprocess import call

print("Warning: This script runs with python3.")
## Check number of input arguments. Should be only 1.
if len(sys.argv) != 2:
    print( 'Number of arguments is ', len(sys.argv)-1, '; one required!.')

# Build folders for processed data
init = initFolderStructure(sys.argv[1])
input_parameters = loadInput(sys.argv[1])

## Time integrate emission spectra
spec = spectra(input_parameters) # (input parameters, i_start, i_end, t_start, t_end) only "input parameters" is mandatory
spec.time_integrate(input_parameters)

## Apply f-scan. Apply either with supergauss parameters or with file path to weights. Make sure spectra are written out before executing this part.
supergauss_parameters = [0.957467, 3.54691, 0.46181, 0.042533, 0.1929, 0.21648] #
spec.fscan(supergauss_parameters)
# spec.fscan(path='/Users/Quincy/Documents/Code/SCFLY_analysis/fscan_weights/fscanweights_Sam')

## Smoothening of spectra
spec.smoothen(np.linspace(1,20,20, dtype='int'), 1.5, 'hanning') # (vector of folders, width [ev], window shape)
spec.smoothen('TOTAL', 1.5, 'hanning') # (vector of folders, width [ev], window shape)

## Extract: temperature-density conditions, population, rates
extract = extract(input_parameters)
extract.temperature_density(input_parameters)
extract.populations(input_parameters,[3,10],'single_ch')
extract.rates(input_parameters,[3,10],'single_ch','coll_ion')

# x=extract.superconfiguration(input_parameters.Z,[3,10],'single_ch')
