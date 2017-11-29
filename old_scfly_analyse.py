#!/usr/bin/env python2.7

#__________________________________________________________________________
#
#  Created by Q.Y. van den Berg on 07/11/2017.
#
#  Process SCFLY output to plot spectra, rates, population, temperature/density
#  python scfly_analyse.py "output directory"
#__________________________________________________________________________

import sys, os, getopt
import numpy
import re
from old_lib_analyse import *
from subprocess import call

# Check number of input arguments. Should be only 1.
if len(sys.argv) != 2:
    print 'Number of arguments is ', len(sys.argv)-1, '; one required!.'

# Build folders for processed data
init = initFolderStructure(sys.argv[1])
input_parameters = loadInput(sys.argv[1])

# Time integrate emission spectra
spec = timeIntSpec(input_parameters) # (input parameters, i_start, i_end, t_start, t_end) only "input parameters" is mandatory
# spec.time_integrate(input_parameters)

# Apply f-scan. Build in synchronisation/lock threads.
# make sure spectra are written out before executing this part.
supergauss_parameters = [0.957467, 3.54691, 0.46181, 0.042533, 0.1929, 0.21648] #
spec.fscan(supergauss_parameters)

# Smoothening of spectra
# spec.broaden(np.linspace(1,20,20), 200, 'GAUSS') # (vector of folders, width [ev], lineshape)
# spec.broaden('TOTAL', 200, 'GAUSS') # (vector of folders, width [ev], lineshape)

# Temperature-density conditions
conditions = extract(input_parameters)
conditions.temperature_density(input_parameters)



    #
    # def __init__(self,filepath):
    #     print os.path.realpath(filepath)+"/processed_data"
    #     print ( os.path.isdir(os.path.realpath(filepath)+"/processed_data")), os.path.dirname(os.path.realpath(filepath))+"/processed_data"
    #
    #     if (not os.path.isdir(os.path.realpath(filepath)+"/processed_data")):
    #         print "Newly created /processed_data directory"
    #         call(["mkdir",os.path.realpath(filepath)+"/processed_data"])
    #     elif (user_yes_no_query("Overwrite old /processed_data folder?")==True):
    #         print "Overwrote old /processed_data directory with new empty directory"
    #         call(["rm",os.path.realpath(filepath)+"/processed_data"])
    #         call(["mkdir",os.path.realpath(filepath)+"/processed_data"])
    #     else:
    #         print "Folder /processed_data already exists and is left untouched."


#
# # Notifiy user of root directory used
# call("clear")
# print "Tree folder set to local directory:"
# call("pwd")
#
# # Load script parameters from file
# input_parameters = loadData(filename)
# input_parameters.displayCount()
#
#
# # Define X-ray beam properties
# lcls_beam = xray_pulse(input_parameters)
#
# # Build directory tree
# masterfolder = input_parameters.material+"/"
# call(["mkdir",masterfolder])
#
# if len(input_parameters.te_range) == 1:
#     energyfolder = str(int(input_parameters.te_range[0]))+"eV_"+str(lcls_beam.cleanIntensity)+"_"+input_parameters.runTAG
# else:
#     energyfolder = str(int(input_parameters.te_range[0]))+"-"+str(int(input_parameters.te_range[-1]))+"eV_"+str(lcls_beam.cleanIntensity)+"_"+input_parameters.runTAG
#
# runfolder = masterfolder+energyfolder
# call(["mkdir",runfolder])
#
#
# # Write intensity grid to file for reference
# T_grid_file = open(runfolder+"/Temperature_grid.dat",'w')
# i = 1
# for temperature in input_parameters.te_range:
#     T_grid_file.write("%1d\t%1.4e\n" %(i, temperature))
#     i += 1
#
# # Build Masterrun file
# mr = open(runfolder+"/masterrun",'w')
# for i in range (0,len(input_parameters.te_range)):
#     mr.write("qsub i"+str(i+1)+"/jobscript_i"+str(i+1)+"\n")
# mr.close()
#
# # Build fixed files
# radiationFile = rad_file(0,input_parameters.Tmax,input_parameters)
#
# runFile = run_file(input_parameters)
#
# for i in range (0,len(input_parameters.te_range)):
#
#     initialFile = ini_file(input_parameters,i)
#
#     folderpath = runfolder+"/i"+str(i+1)+"/"
#     call(["mkdir",folderpath])
#
#     historyFile = hist_file(0,input_parameters.Tmax,input_parameters)
#     historyFile.write(folderpath,input_parameters,i)
#
#     radiationFile.write(folderpath,input_parameters,input_parameters.FEL_photonE,lcls_beam.peakIntensity)
#     initialFile.write(folderpath)
#
#     runFile.write(folderpath,input_parameters,i+1)
#
#     specFile = spec_file()
#     specFile.write_serial(folderpath,input_parameters,i+1)
#     js_file(folderpath,input_parameters,i+1,runfolder)
#
#
# # Save script in right folder for reference
# call("cp "+filename+" "+runfolder,shell=True)
