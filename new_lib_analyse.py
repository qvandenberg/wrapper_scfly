#!/usr/bin/env python2.7

#__________________________________________________________________________
#
#   Library of classes to support scfly_analyse.py. Contains several classes/functions:
#
#   * Time-integrate emission spectra
#   * Extract ion populations vs time
#   * Extract electron temperature and density vs time
#   * F-scan weigh all physical quantities
#   * Plot results
#
#
#__________________________________________________________________________


import sys, os, getopt
import re
import inspect
import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate,arange,array,exp,integrate
from math import *
from tempfile import mkstemp
from shutil import move
from os import remove, close
from subprocess import call
from distutils.util import strtobool
#import matplotlib.pyplot as plt

def user_yes_no_query(question):
    sys.stdout.write('%s [y/n]\n' % question)
    while True:
        try:
            return strtobool(raw_input().lower())
        except ValueError:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')

class initFolderStructure:
    'Creates folders to store processed data.'

    def __init__(self,filepath):
        print os.path.realpath(filepath)+"/processed_data"
        print ( os.path.isdir(os.path.realpath(filepath)+"/processed_data")), os.path.dirname(os.path.realpath(filepath))+"/processed_data"

        if (not os.path.isdir(os.path.realpath(filepath)+"/processed_data")):
            print "Newly created /processed_data directory"
            call(["mkdir",os.path.realpath(filepath)+"/processed_data"])
        elif (user_yes_no_query("Overwrite old /processed_data folder?")==True):
            print "Overwrote old /processed_data directory with new empty directory"
            call(["rm","-r",os.path.realpath(filepath)+"/processed_data"])
            call(["mkdir",os.path.realpath(filepath)+"/processed_data"])
        else:
            print "Folder /processed_data already exists and is left untouched."

# Read input variables. Eventually merge with launcher script.
class loadInput:
    'Common base class for all input parameters'
    fileCount = 0

    def __init__(self,filepath):
        # print "ls: ", os.listdir(os.path.dirname(os.path.abspath(filepath))+"/"+filepath)
        self.basepath = os.path.abspath(filepath)
        # print "List of directories: ",os.listdir(self.basepath)
        # os.listdir(os.path.dirname(os.path.abspath(filepath)))

        for file in os.listdir(self.basepath):
            if re.search('.txt', file):
            #   if len(file)>1:
            #       print "Multiple files with .txt extension present in work folder. Specify input file."
              print "\nInput parameters read from file: ", self.basepath+"/"+file
              with open(os.path.abspath(filepath)+"/"+file,'r') as myfile:
                for line in myfile:
                    li=line.rstrip()
                    if not li.startswith("#"):
                        input_read = line.split()
                        # Read all the parameters from the file

                        if re.match("rootDIR",line):
                            self.rootDIR = input_read[1]
                        elif re.match("SCFLYversion",line):
                            self.SCFLYversion = input_read[1]
                        elif re.match("runTAG",line):
                            self.runTAG = input_read[1]
                        elif re.match("runAUGER",line):
                            self.runAUGER = input_read[1]
                        elif re.match("init_element",line):
                            self.material = input_read[1]
                        elif re.match("init_IPDmodel",line):
                            self.ipd_model = input_read[1]
                        elif re.match("init_ionization",line):
                            self.ionization = float(input_read[1])
                        elif re.match("init_Te",line):
                            self.te = float(input_read[1])
                        elif re.match("init_Ti",line):
                            self.ti = float(input_read[1])
                        elif re.match("init_density",line):
                            self.density = float(input_read[1])
                        elif re.match("fel_energy",line):
                            self.pulse_energy = float(input_read[1])
                        elif re.match("fel_bw",line):
                            self.pulse_bw = float(input_read[1])
                        elif re.match("fel_length",line):
                            self.pulse_length = float(input_read[1])
                        elif re.match("fel_mean",line):
                            self.mean = float(input_read[1])
                        elif re.match("fel_size",line):
                            self.spotsize = float(input_read[1])
                        elif re.match("fel_range",line):
                            self.FEL_photonE = np.linspace(float(input_read[1]),float(input_read[2]),num=int(input_read[3]))
                        elif re.match("init_timesteps",line):
                            self.size_timegrid = int(input_read[1])
                        elif re.match("init_intensitysteps",line):
                            self.size_intensitygrid = int(input_read[1])
                        elif re.match("init_runtime",line):
                            self.Tmax = float(input_read[1])
                        elif re.match("SM_range",line):
                            self.HVrange = np.linspace(float(input_read[1]),float(input_read[2]),num=float(input_read[3]))
                        elif re.match("SM_size",line):
                            self.size = float(input_read[1])

                        loadInput.fileCount += 1

    def displayCount(self):
        print "Total amount of input parameters: %d" % loadInput.fileCount

class timeIntSpec:
    'Routines for time-integration of spectra, f-scan weighing and broadening.'
    # If no time boundaries are set the default is to integrate across all time steps

    def __init__(self,inputParameters, i_start=None,i_end=None,t_start=None,t_end=None):
        # print "input.path %s" % (inputParameters.basepath)
        self.basepath = inputParameters.basepath
        # print 'basepath ', self.basepath
        # folder and log file to store results
        if not os.path.isdir(self.basepath+"/processed_data/spectra"):
            call(["mkdir",self.basepath+"/processed_data/spectra"])
        self.logfile = open(self.basepath+"/processed_data/spectra/intspec.log",'w+a')

        # Determine time and intensity steps to consider, and frequency grid
        self.freq_grid = inputParameters.HVrange

        i_min, i_max = 1, 0
        for directory in os.listdir(inputParameters.basepath):
            if (len(directory)<5): # only act on intensity folders
                i_folder = ((directory.split('i'))) #[-1]
                i_min = min(i_min,int(i_folder[1]))
                i_max = max(i_max,int(i_folder[1]))

        if (i_start == None):
            self.i_start = i_min
        if (i_end == None):
            self.i_end = i_max
        if (t_start == None):
            self.t_start = 2
        if (t_end == None):
            self.t_end = inputParameters.size_timegrid
            self.logfile.write('Analysis run on:\t\t'+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")+"\n")
        self.logfile.write('Intensity folders considered for spectrum: \t i%s until i%s \n' % (str(self.i_start),str(self.i_end)))
        self.logfile.write('Time steps considered for spectrum: \t %s until %s \n' %(str(self.t_start),str(self.t_end)))
        self.logfile.write('Spectral frequency range [eV]: \t %s until %s \n' %(str(min(self.freq_grid)),str(max(self.freq_grid))))
        self.logfile.write('Amount of frequency steps: \t %s \n\n' %(str(len(self.freq_grid))))
        # logfile.close()

    def time_integrate(self, inputParameters): # inputParameters, i_start=None,i_end=None,t_start=None,t_end=None
        # Integrate all intensities independently
        if (os.path.isdir(self.basepath+"/processed_data/spectra")&(user_yes_no_query("Overwrite old /processed_data/spectra folder?")==True)):
            intensity = np.zeros_like(self.freq_grid)
            for i in range(self.i_start, self.i_end+1): # loop over intensity
                i_folderpath = os.path.join(inputParameters.basepath, "i" + str(i))
                data_in_path = os.path.join(i_folderpath, "output")
                # Time integrate
                for j in range(self.t_start, self.t_end+1):
                    data = np.loadtxt(os.path.join(data_in_path, "out.%03d" % (j)), skiprows=2)
                    freqs = data[:,0]
                    idx = 1 # if args.is_intensity else 3
                    intens = data[:,idx]
                    # interpolate the intensity
                    intensity += np.interp(self.freq_grid, freqs, intens)

                # Write out intensities to file

                file_specout = open(self.basepath+"/processed_data/spectra/time-integrated_i"+str(i),'w+a')
                file_specout.write("%s\t%s\n" %('E [eV]','Intensity'))
                for k in range(len(self.freq_grid)):
                    file_specout.write("%1.2f %1.4e\n" %(self.freq_grid[k],intensity[k]))
                    # print ("%1.2f %1.4e\n" %(self.freq_grid[k],intensity[k]))
                file_specout.close()
            print "Time-integrated spectra written to folder ",self.basepath+"/processed_data/spectra/"
        else:
            print "Time-integrated spectra already exist in folder %s." %(self.basepath+"/processed_data/spectra")

            # Plot for inspection
            # plt.plot(self.freq_grid, intensity)
            # plt.show()

    def fscan(self, supergauss_parameters, i_start = None, i_end=None):
# generalise to do the same for temperature, density, population
        if i_start == None:
            i_start = self.i_start
        if i_end == None:
            i_end = self.i_end
        ## Spectrum
        intensity_grid = np.loadtxt(os.path.join(self.basepath, 'intensity_grid.dat'), skiprows=0)
        intensity_grid[:,0]= [int(i) for i in intensity_grid[:,0]] # map(int,intensity_grid[:,0])
        intensity_grid[:,1]= intensity_grid[:,1]/max(intensity_grid[:,1])

        # Calculate super-gauss function on a grid
        S_fscan = np.linspace(0,500,25000)
        I_fscan_unnormalised = supergauss_parameters[0]*np.exp(-np.power(S_fscan/supergauss_parameters[1],supergauss_parameters[2]))
        + supergauss_parameters[3]*np.exp(-np.power(S_fscan/supergauss_parameters[4],supergauss_parameters[5]))
        I_fscan = I_fscan_unnormalised/max(I_fscan_unnormalised)

        # interpolate simulated intensity_grid onto super-gauss function to obtain weights
        S_interpolated = np.interp(intensity_grid[:,1], I_fscan[::-1], S_fscan[::-1]) # swap direction x, y to make them increasing
        surface_weights = np.convolve(S_interpolated, np.ones((2,))/2)[(2-1):]
        surface_weights = np.diff(S_interpolated,n=1) # this operation reduces array length by one
        surface_weights = np.append(surface_weights, 2*surface_weights[-1]-surface_weights[-2]) # duplicate last item to restore length with lin extrapolation
        surface_weights = surface_weights/sum(surface_weights)

        # weigh spectra to form total spectrum
        total_spectrum = np.zeros((len(self.freq_grid),))
        for i in range(len(surface_weights)):
            datapath = self.basepath+"/processed_data/spectra/time-integrated_i"+str(int(i+1))
            data = np.loadtxt(datapath, skiprows=1)
            x = data[:,0] # energy, eV
            y = data[:,1] # intensity
            # Plot for inspection
            # if i==10:
            # plt.plot(x, y)
            # plt.show()
            total_spectrum += surface_weights[i]*np.interp(self.freq_grid,x,y)

        # Write weights to file
        file_weights = open(self.basepath+"/processed_data/fscanweights",'w+a')
        file_weights.write("%s\t%s\t%s\n" %('Index','Intensity','Weight'))
        for k in range(len(surface_weights)):
            file_weights.write("%d %1.4e %1.3f\n" %(k+1, I_fscan_unnormalised[k], surface_weights[k]))

        # Write spectrum to file
        file_specout = open(self.basepath+"/processed_data/spectra/time-integrated_total",'w+a')
        file_specout.write("%s\t%s\n" %('E [eV]','Intensity'))
        for k in range(len(self.freq_grid)):
            file_specout.write("%1.2f %1.4e\n" %(self.freq_grid[k],total_spectrum[k]))
        file_specout.close()

        # Update log file
        self.logfile.write('Constructed total spectrum from f-scan weighted spectra.\n')
        self.logfile.write("f-scan super-gauss parameters used are [%1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f] \n" %(supergauss_parameters[0], supergauss_parameters[1],
        supergauss_parameters[2], supergauss_parameters[3], supergauss_parameters[4], supergauss_parameters[5]))
        self.logfile.write("f-scan weights used:\n")
        self.logfile.write('\t'.join(map(str,np.around(surface_weights,5))))
        self.logfile.write("\n")

        ## Plot for inspection
        # plt.plot(self.freq_grid, total_spectrum)
        # plt.show()

    def broaden(self, i_folders, sigma, line_shape):
        # Somehow this results in negative spectra. Fix first before applying.
        print "Fix broaden/smoothen bug first before using further."
        if (isinstance(i_folders,(list, tuple, np.ndarray))==True):
            for i in i_folders:
                # print 'folders', i
                # Save copy of original spectrum and update log file
                if not os.path.exists(self.basepath+"/processed_data/spectra/oldcopies"):
                    print 'Created directory for old copies of spectra'
                    call(["mkdir",self.basepath+"/processed_data/spectra/oldcopies"])
                # print 'line 179 ',self.basepath+"/processed_data/spectra/oldcopies/time-integrated_i"+str(int(i))+'_'+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S").replace(' ','_')
                call("cp "+self.basepath+"/processed_data/spectra/time-integrated_i"+str(int(i)) +' '+self.basepath+"/processed_data/spectra/oldcopies/time-integrated_i"+str(int(i))+'_'+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S").replace(' ','_'),shell=True)
                self.logfile.write('Broadened i'+str(int(i))+' spectrum by a %s line shape with width: %1.2f \n' %(line_shape, sigma))
                self.logfile.write('Copy of i'+str(int(i))+' spectrum made with timestamp '+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")+ "\n" %())

                # Read spectra from file
                datapath = self.basepath+"/processed_data/spectra/time-integrated_i"+str(int(i))
                data = np.loadtxt(datapath, skiprows=1)
                x = data[:,0] # energy, eV
                y = data[:,1] # intensity

                # Construct line shape
                L = max(x)-min(x)
                x_convolve = np.linspace(-L/2,L/2,len(x))

                if line_shape == 'GAUSS':
                    y_convolve = lambda x: np.exp(-x**2/(2*sigma**2))
                    norm = integrate.quad(y_convolve,-L/2,L/2)[0]
                else:
                    print "Currently only supports Gaussian lineshape. Supply argument \'GAUSS\'"
                # Convolve original spectrum with line shape for smoothening
                y_convolve = np.exp(-x_convolve**2/(2*sigma**2))/norm
                y_smooth = np.convolve(y, y_convolve, mode='same')
                y_smooth = y_smooth*np.trapz(x,y)/np.trapz(x,y_smooth)
                # Write result to file
                file_specout = open(datapath,'w+a')
                file_specout.write("%s\t%s\n" %('E [eV]','Intensity'))
                for k in range(len(x)):
                    file_specout.write("%1.2f %1.2f\n" %(x[k],y_smooth[k]))
                file_specout.close()

        elif (i_folders == 'TOTAL'):
            print 'total spectrum'
            # Save copy of original spectrum and update log file
            if not os.path.exists(self.basepath+"/processed_data/spectra/oldcopies"):
                print 'Created directory for old copies of spectra'
                call(["mkdir",self.basepath+"/processed_data/spectra/oldcopies"])
            call("cp "+self.basepath+"/processed_data/spectra/time-integrated_total "+' '+self.basepath+"/processed_data/spectra/oldcopies/time-integrated_total_"+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S").replace(' ','_'),shell=True)
            self.logfile.write('Broadened total spectrum by a %s line shape with width: %1.2f \n' %(line_shape, sigma))
            self.logfile.write('Copy of original spectrum made with timestamp '+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")+ "\n")

            # Read spectra from file
            datapath = self.basepath+"/processed_data/spectra/time-integrated_total"
            data = np.loadtxt(datapath, skiprows=1)
            x = data[:,0] # energy, eV
            y = data[:,1] # intensity

            # Construct line shape
            L = max(x)-min(x)
            x_convolve = np.linspace(-L/2,L/2,len(x))

            if line_shape == 'GAUSS':
                self.line_shape = 'GAUSS'
                y_convolve = lambda x: np.exp(-x**2/(2*sigma**2))
                norm = integrate.quad(y_convolve,-L/2,L/2)[0]
            else:
                print "Currently only supports Gaussian lineshape. Supply argument \'GAUSS\'"
            y_convolve = np.exp(-x_convolve**2/(2*sigma**2))/norm
            y_smooth = np.convolve(y, y_convolve, mode='same')
            y_smooth = y_smooth*np.trapz(x,y)/np.trapz(x,y_smooth)
            # Write result to file
            file_specout = open(datapath,'w+a')
            file_specout.write("%s\t%s\n" %('E [eV]','Intensity'))
            for k in range(len(x)):
                file_specout.write("%1.2f %1.2f\n" %(x[k],y_smooth[k]))
            file_specout.close()

        else:
            print "Either operate on a range of input directories, or the total f-scan weighted spectrum. Check input arguments."

class extract:
    'Extract populations, temperature, density'

    def __init__(self,inputParameters, i_start=None,i_end=None,t_start=None,t_end=None):
        self.basepath = inputParameters.basepath

        # Determine time and intensity steps to consider
        i_min, i_max = 1, 0
        for directory in os.listdir(inputParameters.basepath):
            if (len(directory)<5): # only act on intensity folders
                i_folder = ((directory.split('i'))) #[-1]
                i_min = min(i_min,int(i_folder[1]))
                i_max = max(i_max,int(i_folder[1]))
        if (i_start == None):
            self.i_start = i_min
        if (i_end == None):
            self.i_end = i_max
        if (t_start == None):
            self.t_start = 2
        if (t_end == None):
            self.t_end = inputParameters.size_timegrid

        self.conditions_logfile = open(self.basepath+"/processed_data/conditions/conditions.log",'w+a')
        if ((not os.path.isdir(self.basepath+"/processed_data/conditions")) or (user_yes_no_query("Overwrite /processed_data/conditions extracted plasma conditions folder?") == True)):
            call(["rm","-r",os.path.join(self.basepath,"processed_data/conditions")])
            call(["mkdir",self.basepath+"/processed_data/conditions"])
            print "New folder created: %s" %(os.path.join(self.basepath,"processed_data/conditions"))
        else:
            print "Folder already existed to store conditions: %s" %(os.path.join(self.basepath,"processed_data/conditions"))

    def temperature_density(self,inputParameters):
        self.conditions_logfile.write('Intensity folders considered for temperature-density evolution:\ti%s until i%s\n' % (str(self.i_start),str(self.i_end)))
        self.conditions_logfile.write('Time steps considered:\t%s until %s\n' %(str(self.t_start),str(self.t_end)))

        for i in range(self.i_start,self.i_end+1):
            i_folderpath = os.path.join(inputParameters.basepath, "i" + str(i),"output","zb.i"+str(i))
            data = np.loadtxt(i_folderpath, skiprows=1)
            trho_out = open(self.basepath+"/processed_data/conditions/trho_i%d.txt" %(i),'w+a')
            trho_out.write("%s\t%s\t%s\t%s\n" %('Time index','Time [s]','Temperature [eV]','Density [/cc]'))
            # Write out data
            for j in range(data.shape[0]):
                trho_out.write("%d\t%1.2e\t%1.2f\t%1.2e\n" %(j,data[j,1],data[j,2],data[j,3])) # index, time, temperature, density

        # Weigh temperatures together to effective temperature with fscan
        weights = []
        if os.path.isfile(self.basepath+"/processed_data/fscanweights")==True:
            data = np.loadtxt(self.basepath+"/processed_data/fscanweights", skiprows=1)
            weights = data[:,2]
        trhofscan_out = open(self.basepath+"/processed_data/conditions/trho_fscanweighted.txt",'w+a')
        trhofscan_out.write("%s\t%s\t%s\t%s\n" %('Time index','Time [s]','Temperature [eV]','Density [/cc]'))

        time_fscan =
        T_fscan = time_fscan
        rho_fscan = time_fscan

        for i in range(self.i_start,self.i_end+1):
            data = np.loadtxt(self.basepath+"/processed_data/conditions/trho_i%d.txt" %(i), skiprows=1)









class xray_pulse:
    'Defines the X-ray intensities needed for Scfly'

    def __init__(self,parameters):
        self.peakIntensity = parameters.pulse_energy/(parameters.pulse_length*parameters.spotsize)

        power = floor(log10(self.peakIntensity))
        coeff =  round(self.peakIntensity/10**power)
        self.cleanIntensity = coeff*10**power


    # This function generates a good set of calculations which can then be fit
    # to a generic spot_size grid.
    # Extends down from peak intensity 5 orders of magnitude if 20 steps used.
    def init_Intgrid (self,parameters,length,sigma):
        xlist = np.linspace(0.01,length,num=parameters.size_intensitygrid)
        LogGrid = [exp(-(xvar/sigma)) for xvar in xlist]
        self.intensityGrid = [self.peakIntensity * gridpoints for gridpoints in LogGrid]

    def init_Intgrid_Gauss (self,parameters,length,sigma):
        xlist = np.linspace(0.01,length,num=parameters.size_intensitygrid)
        gaussGrid = [exp(-0.5 * (xvar/sigma)**2) for xvar in xlist]
        self.intensityGrid = [self.peakIntensity * gridpoints for gridpoints in gaussGrid]



class hist_file:
    'Writes the history file for Scfly'

    def __init__(self,start_time,end_time,parameters):
        self.time_vector = 1e-15 * np.linspace(start_time,end_time,num=parameters.size_timegrid)

    def write(self,folder,parameters):
        hf = open(folder+"/hist",'w')
        hf.write("time rho size te\n")

        for i in range (0,parameters.size_timegrid):
            hf.write("%1.4e %1.4e %1.4e %1.4e \n" %(self.time_vector[i],parameters.density,parameters.size*1e-4,parameters.te))

        hf.close()


class rad_file:
    'Writes the radiation file for Scfly'

    def __init__(self,start_time,end_time,parameters):
        #step = (end_time - start_time) / (parameters.size_timegrid-1)
        self.time_vector = 1e-15 * np.linspace(start_time,end_time,num=parameters.size_timegrid)
        self.sigma = parameters.pulse_length/(2*sqrt(2*log(2)))
        try:
            self.mean = parameters.mean
        except NameError:
            self.mean =  1e-15 * (start_time + end_time)/2



    def write(self,folder,parameters,photon_energy,intensity):
        # Beam parameters in converted units
        bandwidth = parameters.pulse_bw * 0.01 * photon_energy;     # in eV
        bw = bandwidth * 2.418e14;                                  # in Hz
        erg = 1e-7;                                                 # in Joules
        Peak_Fluence = intensity / (erg * bw * 4*pi)                # in scfly units

        ## For a Gaussian temporal profile
        FEL_timeGrid = [Peak_Fluence * exp(-0.5 * ((tvar-self.mean)/self.sigma)**2) for tvar in self.time_vector]

        sigma_bw = bandwidth/(2*sqrt(2*log(2)));
        energy_vector = np.linspace(photon_energy-5*sigma_bw,photon_energy+5*sigma_bw,num=11)
        FEL_bandwidth = [exp(-0.5 * ((bvar-photon_energy)/sigma_bw)**2) for bvar in energy_vector]

        rf = open(folder+"/hv_file",'w')
        rf.write("%d %d\n" %(parameters.size_timegrid,len(energy_vector)))

        for i in range (0,len(energy_vector)):
            rf.write("%1.2f " %(energy_vector[i]))
        rf.write("\n")

        for i in range (0,len(self.time_vector)):
            rf.write("%1.3e 0.000e+00\n" %(self.time_vector[i]))
            for j in range (0,len(FEL_bandwidth)):
                rf.write("%1.3e " %(FEL_timeGrid[i] * FEL_bandwidth[j]))
            rf.write("\n")
        rf.close()


class ini_file:
    'Writes the initial file for Scfly'

    def __init__(self,parameters):
        # This version only supports C, Al, Mg, Fe and Cu with integer charge state
        self.mixture = 'no';

        if parameters.material == "C":

            self.Z = 6;
            # Force charge state to be integer
            charge = 6 - round(parameters.ionization)

            if charge == 1:
                self.atom_name = 'h_100001';
            elif charge == 2:
                self.atom_name = 'he200002';
            elif charge == 3:
                self.atom_name = 'li210002';
            elif charge == 4:
                self.atom_name = 'be220002';
            elif charge == 5:
                self.atom_name = 'b_230002';
            elif charge == 6:
                self.atom_name = 'c_240002';
            else:
                print "Warning: Incorrect charge state.\nPlease choose a charge state between 0 and 5 for Carbon!"


        if parameters.material == "Mg":

            self.Z = 12;
            # Force charge state to be integer
            charge = 12 - round(parameters.ionization)

            if charge == 1:
                self.atom_name = 'h_100001';
            elif charge == 2:
                self.atom_name = 'he200002';
            elif charge == 3:
                self.atom_name = 'li210002';
            elif charge == 4:
                self.atom_name = 'be220002';
            elif charge == 5:
                self.atom_name = 'b_230002';
            elif charge == 6:
                self.atom_name = 'c_240002';
            elif charge == 7:
                self.atom_name = 'n_250002';
            elif charge == 8:
                self.atom_name = 'o_260002';
            elif charge == 9:
                self.atom_name = 'f_270002';
            elif charge == 10:
                self.atom_name = 'ne280002';
            elif charge == 11:
                self.atom_name = 'na281003';
            elif charge == 12:
                self.atom_name = 'mg282003';
            else:
                print "Warning: Incorrect charge state.\nPlease choose a charge state between 0 and 11 for Mg!"


        if parameters.material == "MgF2":

            self.mixture = 'yes';
            # This bit has not been properly coded yet - only hardcoded at the bottom for now.

        if parameters.material == "Al":

            self.Z = 13;
            # Force charge state to be integer
            charge = 13 - round(parameters.ionization)

            if charge == 1:
                self.atom_name = 'h_100001';
            elif charge == 2:
                self.atom_name = 'he200002';
            elif charge == 3:
                self.atom_name = 'li210002';
            elif charge == 4:
                self.atom_name = 'be220002';
            elif charge == 5:
                self.atom_name = 'b_230002';
            elif charge == 6:
                self.atom_name = 'c_240002';
            elif charge == 7:
                self.atom_name = 'n_250002';
            elif charge == 8:
                self.atom_name = 'o_260002';
            elif charge == 9:
                self.atom_name = 'f_270002';
            elif charge == 10:
                self.atom_name = 'ne280002';
            elif charge == 11:
                self.atom_name = 'na281003';
            elif charge == 12:
                self.atom_name = 'mg282003';
            elif charge == 13:
                self.atom_name = 'al283003';

            else:
                print "Warning: Incorrect charge state.\nPlease choose a charge state between 0 and 12 for Al!"


        if parameters.material == "Fe":

            self.Z = 26;
            # Force charge state to be integer
            charge = 26 - round(parameters.ionization)

            if charge == 1:
                self.atom_name = 'h_100001';
            elif charge == 2:
                self.atom_name = 'he200002';
            elif charge == 3:
                self.atom_name = 'li210002';
            elif charge == 4:
                self.atom_name = 'be220002';
            elif charge == 5:
                self.atom_name = 'b_230002';
            elif charge == 6:
                self.atom_name = 'c_240002';
            elif charge == 7:
                self.atom_name = 'n_250002';
            elif charge == 8:
                self.atom_name = 'o_260002';
            elif charge == 9:
                self.atom_name = 'f_270002';
            elif charge == 10:
                self.atom_name = 'ne280002';
            elif charge == 11:
                self.atom_name = 'na281003';
            elif charge == 12:
                self.atom_name = 'mg282003';
            elif charge == 13:
                self.atom_name = 'al283003';
            elif charge == 14:
                self.atom_name = 'si284003';
            elif charge == 15:
                self.atom_name = 'p_285003';
            elif charge == 16:
                self.atom_name = 's_286003';
            elif charge == 17:
                self.atom_name = 'cl287003';
            elif charge == 18:
                self.atom_name = 'ar288003';
            elif charge == 19:
                self.atom_name = 'k_289003';
            elif charge == 20:
                self.atom_name = 'ca28a003';
            elif charge == 21:
                self.atom_name = 'sc28b003';
            elif charge == 22:
                self.atom_name = 'ti28c003';
            elif charge == 23:
                self.atom_name = 'v_28d003';
            elif charge == 24:
                self.atom_name = 'cr28e003';
            elif charge == 25:
                self.atom_name = 'mn28f003';

            else:
                print "Warning: Incorrect charge state.\nPlease choose a charge state between 4 and 25 for Fe!"


        if parameters.material == "Cu":

            self.Z = 29;
            # Force charge state to be integer
            charge = 29 - round(parameters.ionization)

            if charge == 1:
                self.atom_name = 'h_100001';
            elif charge == 2:
                self.atom_name = 'he200002';
            elif charge == 3:
                self.atom_name = 'li210002';
            elif charge == 4:
                self.atom_name = 'be220002';
            elif charge == 5:
                self.atom_name = 'b_230002';
            elif charge == 6:
                self.atom_name = 'c_240002';
            elif charge == 7:
                self.atom_name = 'n_250002';
            elif charge == 8:
                self.atom_name = 'o_260002';
            elif charge == 9:
                self.atom_name = 'f_270002';
            elif charge == 10:
                self.atom_name = 'ne280002';
            elif charge == 11:
                self.atom_name = 'na281003';
            elif charge == 12:
                self.atom_name = 'mg282003';
            elif charge == 13:
                self.atom_name = 'al283003';
            elif charge == 14:
                self.atom_name = 'si284003';
            elif charge == 15:
                self.atom_name = 'p_285003';
            elif charge == 16:
                self.atom_name = 's_286003';
            elif charge == 17:
                self.atom_name = 'cl287003';
            elif charge == 18:
                self.atom_name = 'ar288003';
            elif charge == 19:
                self.atom_name = 'k_289003';
            elif charge == 20:
                self.atom_name = 'ca28a003';
            elif charge == 21:
                self.atom_name = 'sc28b003';
            elif charge == 22:
                self.atom_name = 'ti28c003';
            elif charge == 23:
                self.atom_name = 'v_28d003';
            elif charge == 24:
                self.atom_name = 'cr28e003';
            elif charge == 25:
                self.atom_name = 'mn28f003';
            elif charge == 26:
                self.atom_name = 'fe28g003';
            elif charge == 27:
                self.atom_name = 'co28h003';
            elif charge == 28:
                self.atom_name = 'ni28i003';

            else:
                print "Warning: Incorrect charge state.\nPlease choose a charge state between 2 and 28 for Cu!"

    def write(self,folder):
        intf = open(folder+"/initial.dat",'w')
        if self.mixture == 'no':
            intf.write(self.atom_name+" 1.d0\n")
        else:
            intf.write("12\tne280002\t1.d0\n")
            intf.write("9\tf_270002\t2.d0\n")
        intf.close()

class spec_file:
    'Writes spectral model file for Scfly'

#    def __init__(self,folder,parameters,index):
#        spf = open(folder+"/runspect",'w')
#        spf.write("i"+str(index)+"\n")
#        spf.write("2 %d 1\n" %(parameters.size_timegrid+1))
#        spf.write("%1.1f %1.1f\n" %(parameters.HVrange[0],parameters.HVrange[-1]))
#        spf.write("%1.2e\n" %(parameters.size*1e-4))
#        spf.write("0.000\n2.000\n")
#        spf.close()

    def write_serial(self,folder,parameters,index):

        if parameters.material == "C":
            self.Z = 6;
        elif parameters.material == "Mg":
            self.Z = 12;
        elif parameters.material == "Al":
            self.Z = 13;
        elif parameters.material == "Fe":
            self.Z = 26;
        elif parameters.material == "Cu":
            self.Z = 29;
        else:
            print "Warning: Element is not supported! Please use C, Mg, Al, Fe or Cu."

        spf = open(folder+"/runspect",'w')
        spf.write("scfly-output-file:\t\t"+str(self.Z)+".i"+str(index)+"\n")
        spf.write("ipd-file:\t\t\t"+"i"+str(index)+".ipd\n")
        spf.write("transition-binary-file:\t\t"+str(self.Z)+".trb\n")
        spf.write("output-directory:\t\t\toutput/\n")
        spf.write("output-type\t\ttext\n")
        spf.write("time-grid:\t\t\t2 %d 1\n" %(parameters.size_timegrid))
        spf.write("energy-range:\t\t %1.1f %1.1f\n" %(parameters.HVrange[0],parameters.HVrange[-1]))
        spf.write("size:\t\t\t\t%1.2e\n" %(parameters.size*1e-4))
        spf.write("fwhm-gaussian:\t\t\t%1.2e\n" %(0.0))
        spf.write("fwhm-lorentzian:\t\t%1.2e\n" %(0.0))
        spf.close()


class run_file:
    'Writes runfile for Scfly'

    def __init__(self,parameters):
        if parameters.material == "C":
            self.Z = 6;
        elif parameters.material == "Mg":
            self.Z = 12;
        elif parameters.material == "Al":
            self.Z = 13;
        elif parameters.material == "Fe":
            self.Z = 26;
        elif parameters.material == "Cu":
            self.Z = 29;
        else:
            print "Warning: Element is not supported! Please use C, Mg, Al, Fe or Cu."

    def write(self,folder,parameters,index):
        runf = open(folder+"/runfile_i"+str(index),'w')
        runf.write("z "+str(self.Z)+"\n")
        runf.write("tr trfile hv_file\n")
        runf.write("ti file\n")
        runf.write("opacity file\n")
        runf.write("evolve td\n")
        runf.write("initial file initial.dat\n")
        runf.write("history hist rho\n")
        runf.write("io 2\n") # to print all collision rates

        #runf.write("lattice 1.0\n")
        if parameters.ipd_model == "EK":
            runf.write("ipd ek\n")
        elif parameters.ipd_model == "SP":
            runf.write("ipd sp\n")
        elif parameters.ipd_model == "experimental":
            runf.write("ipd file\n")
        elif parameters.ipd_model == "atomic":
            runf.write("ipd file\n")
        else:
            print "Warning: IPD model not supported. Please specify EK, SP, experimental or atomic."

        runf.write("outfile i"+str(index)+"\n")
        runf.write("isos 1 "+str(self.Z)+" 4\n")
        runf.write("end\n")
        runf.close()

    def writeAuger(self,folder,parameters,index):
        runf = open(folder+"/runfile_i"+str(index),'w')
        runf.write("z "+str(self.Z)+"\n")
        runf.write("tr trfile hv_file\n")
        runf.write("ti file\n")
        runf.write("opacity file\n")
        runf.write("evolve td\n")
        runf.write("initial file initial.dat\n")
        runf.write("fe auge\n")
        runf.write("history hist rho\n")
        runf.write("io 2\n") # to print all collision rates

        #runf.write("lattice 1.0\n")
        if parameters.ipd_model == "EK":
            runf.write("ipd ek\n")
        elif parameters.ipd_model == "SP":
            runf.write("ipd sp\n")
        else:
            print "Warning: IPD model not supported. Please specify EK or SP."

        runf.write("outfile i"+str(index)+"\n")
        runf.write("isos 1 "+str(self.Z)+" 10\n")
        runf.write("end\n")
        runf.close()



class js_file:
    'Write the jobscript file for running scfly and spec within the intensity folder ixx'

    def __init__(self,folder,parameters,index,basefolder):
        jsf = open(folder+"/jobscript_i"+str(index)+".sh",'w')
        jsf.write("#!/bin/bash\n\n")
        jsf.write("PERMDIR="+parameters.rootDIR+basefolder+"/i"+str(index)+"\n")
        #jsf.write("\nEXECDIR="+parameters.rootDIR+"scfly_exec/"+parameters.ipd_model+"\n")
        jsf.write("ATOMDIR="+parameters.rootDIR+"scfly_executables/ATOM\n")
        jsf.write("EXECDIR="+parameters.rootDIR+"scfly_executables/EXEC\n")
#        jsf.write("SPECDIR="+parameters.rootDIR+"scfly_executables/SPEC\n")

        if parameters.material == "C":
            jsf.write("SPECDIR="+parameters.rootDIR+"scfly_executables/spectral/files/transition/06-C\n\n")
            self.Z = 6;
        elif parameters.material == "Mg":
            jsf.write("SPECDIR="+parameters.rootDIR+"scfly_executables/spectral/files/transition/12-Mg\n\n")
            self.Z = 12;
        elif parameters.material == "Al":
            jsf.write("SPECDIR="+parameters.rootDIR+"scfly_executables/spectral/files/transition/13-Al\n\n")
            self.Z = 13;
        elif parameters.material == "Fe":
            jsf.write("SPECDIR="+parameters.rootDIR+"scfly_executables/spectral/files/transition/26-Fe\n\n")
            self.Z = 26;
        elif parameters.material == "Cu":
            jsf.write("SPECDIR="+parameters.rootDIR+"scfly_executables/spectral/files/transition/29-Cu\n\n")
            self.Z = 29;
        else:
            print "Warning: Element is not supported! Please use C, Mg, Al, Fe or Cu."

        jsf.write("cd $EXECDIR/Version_"+parameters.SCFLYversion+"\n")
        jsf.write("make\n\n")
        jsf.write("cd $PERMDIR\n\n")
        jsf.write("cp $EXECDIR/Version_"+parameters.SCFLYversion+"/scfly $PERMDIR\n")
#        jsf.write("cp $EXECDIR/Spectral_module/spec $PERMDIR\n")
#        jsf.write("cp $SPECDIR/* $PERMDIR\n")
        jsf.write("cp $ATOMDIR/"+str(self.Z)+".ipd $PERMDIR\n")
        jsf.write("cp $ATOMDIR/atomic.data $PERMDIR\n")
        jsf.write("cp $ATOMDIR/atomic.inp."+str(self.Z)+" $PERMDIR\n\n")
        jsf.write("./scfly runfile_i"+str(index)+"\n\n")

        jsf.write("mkdir $PERMDIR/"+str(self.Z)+"\n")
        if parameters.ipd_model == "experimental":
            if (len(str(self.Z))==1):
                jsf.write("cp $ATOMDIR/0"+str(self.Z)+"fixed.ipd"+" $PERMDIR/"+str(self.Z)+"/i"+str(index)+".ipd\n")
            elif (len(str(self.Z))==2):
                 jsf.write("cp $ATOMDIR/"+str(self.Z)+"fixed.ipd"+" $PERMDIR/"+str(self.Z)+"/i"+str(index)+".ipd\n")
            else:
                 print "Warning: Experimental IPD file not available."

        jsf.write("cp $PERMDIR/"+str(self.Z)+".i"+str(index)+" $PERMDIR/"+str(self.Z)+"\n")
        jsf.write("cp $PERMDIR/i"+str(index)+".ipd $PERMDIR/"+str(self.Z)+"\n")
        jsf.write("cp $SPECDIR/"+str(self.Z)+".trb $PERMDIR/"+str(self.Z)+"\n")
        jsf.write("cp $EXECDIR/spec $PERMDIR/"+str(self.Z)+"\n")
        jsf.write("cp $PERMDIR/runspect $PERMDIR/"+str(self.Z)+"\n")
        jsf.write("cd "+str(self.Z)+"\n")
        jsf.write("./spec runspect\n")
        jsf.write("cd ..\n")

        jsf.write("mkdir $PERMDIR/output\n")
        jsf.write("mv $PERMDIR/"+str(self.Z)+"/"+str(self.Z)+".i"+str(index)+" $PERMDIR/output\n")
        jsf.write("mv $PERMDIR/"+str(self.Z)+"/i"+str(index)+".ipd $PERMDIR/output\n")
        jsf.write("mv $PERMDIR/"+str(self.Z)+"/output/out.* $PERMDIR/output\n")
        jsf.write("mv $PERMDIR/rate.* $PERMDIR/output\n")
        jsf.write("mv $PERMDIR/zb* $PERMDIR/output\n")
        jsf.write("mv $PERMDIR/bb* $PERMDIR/output\n")
        jsf.write("mv $PERMDIR/bf* $PERMDIR/output\n")

        jsf.write("rm ./*.trb\nrm ./spec\n")
        jsf.write("rm ./*.o\nrm ./scfly\nrm ./atomic.data\n")
        jsf.write("rm ./atomic.inp.*\n")
        jsf.write("rm ./calcu.dat\n")
        jsf.write("rm -r $PERMDIR/"+str(self.Z) +"\n")
        jsf.write("rm ./12.0*\n")
        jsf.write("rm ./0*.ipd\n")
        jsf.write("rm ./1*.ipd\n")
        jsf.write("rm ./readme\n\n")

        jsf.write("tar -cjvf spec.tbz  sp*\n")
        jsf.close()


class remote_scfly:

    def __init__ (self, filename):

        # Load script parameters from file
        parameters = loadData(filename)

        # Define X-ray beam properties
        lcls_beam = xray_pulse(parameters)
        lcls_beam.init_Intgrid(parameters,20,4)
        self.cleanINT = lcls_beam.cleanIntensity
        self.intensityGrid = lcls_beam.intensityGrid

        # Build directory tree
        self.masterfolder = parameters.material+"/"

        if len(parameters.FEL_photonE) == 1:
            self.energyfolder = str(int(parameters.FEL_photonE[0]))+"eV_"+str(lcls_beam.cleanIntensity)+"_"+parameters.runTAG
        else:
            self.energyfolder = str(int(parameters.FEL_photonE[0]))+"-"+str(int(parameters.FEL_photonE[-1]))+"eV_"+str(lcls_beam.cleanIntensity)+"_"+parameters.runTAG

        self.runfolder = self.masterfolder+self.energyfolder
        self.rootDIR  = parameters.rootDIR
        self.FEL_grid = parameters.FEL_photonE
        self.INT_grid = parameters.size_intensitygrid
        self.TIM_grid = parameters.size_timegrid
        self.wl_range = parameters.HVrange
        self.material = parameters.material

    def downloadBase(self):
        print 'Downloading key data from cplxint1.physics.ox.ac.uk...'

        if self.material == "C":
            self.Z = 6;
        elif self.material == "Mg":
            self.Z = 12;
        elif self.material == "Al":
            self.Z = 13;
        elif self.material == "Fe":
            self.Z = 26;
        elif self.material == "Cu":
            self.Z = 29;
        else:
            print "Warning: Element is not supported! Please use C, Mg, Al, Fe or Cu."

        for i in range (0,len(self.FEL_grid)):
            for j in range (0,self.INT_grid):

                folderpath = self.rootDIR+self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)+"/output/"
                destination_folder = self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)+"/output/"
                #call('mkdir '+destination_folder,shell=True)
                call('rsync --progress -avz --rsh=ssh -r vandenberg@cplxint1.physics.ox.ac.uk:'+folderpath+str(self.Z)+".i"+str(i*int(self.INT_grid)+j+1)+" "+destination_folder,shell=True)
                call('rsync --progress -avz --rsh=ssh -r vandenberg@cplxint1.physics.ox.ac.uk:'+folderpath+"zb.i"+str(i*int(self.INT_grid)+j+1)+" "+destination_folder,shell=True)
                #call('rsync --progress -avz --rsh=ssh -r vandenberg@cplxint1.physics.ox.ac.uk:'+folderpath+"pp.i"+str(i*int(self.INT_grid)+j+1)+" "+destination_folder,shell=True)
                #call('rsync --progress -avz --rsh=ssh -r vandenberg@cplxint1.physics.ox.ac.uk:'+folderpath+"popout.d "+destination_folder,shell=True)
                call('rsync --progress -avz --rsh=ssh -r vandenberg@cplxint1.physics.ox.ac.uk:'+folderpath+"i"+str(i*int(self.INT_grid)+j+1)+".ipd "+destination_folder,shell=True)
                #call('rsync --progress -avz --rsh=ssh -r vandenberg@cplxint1.physics.ox.ac.uk:'+folderpath+"spec.tbz "+destination_folder,shell=True)

                # Unpack spectral tarball and remove
                call('tar xjf '+destination_folder+'spec.tbz -C '+destination_folder,shell=True)

        print 'Download completed.\n'

    def downloadOpacity(self):
        print 'Downloading opacity data from cplxint1.physics.ox.ac.uk...'
        for i in range (0,len(self.FEL_grid)):
            for j in range (0,self.INT_grid):

                folderpath = self.rootDIR+self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)+"/output/"
                destination_folder = self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)+"/output/"
                #call('mkdir '+destination_folder,shell=True)
                call('rsync --progress -avz --rsh=ssh -r vandenberg@cplxint1.physics.ox.ac.uk:'+folderpath+"op* "+destination_folder,shell=True)

        print 'Download completed.\n'

    def downloadRates(self):
        print 'Downloading rate files from cplxint1.physics.ox.ac.uk...'
        for i in range (0,len(self.FEL_grid)):
            for j in range (0,self.INT_grid):

                folderpath = self.rootDIR+self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)+"/output/"
                destination_folder = self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)+"/output/"
                #call('mkdir '+destination_folder,shell=True)
                call('rsync --progress -avz --rsh=ssh -r vandenberg@cplxint1.physics.ox.ac.uk:'+folderpath+"rate.* "+destination_folder,shell=True)
                call('rsync --progress -avz --rsh=ssh -r vandenberg@cplxint1.physics.ox.ac.uk:'+folderpath+"sp0* "+destination_folder,shell=True)
                call('rsync --progress -avz --rsh=ssh -r vandenberg@cplxint1.physics.ox.ac.uk:'+folderpath+"op0* "+destination_folder,shell=True)

        print 'Download of rate files completed.\n'


    def downloadFulltree(self):
        print 'Warning: Downloading full data tree from cplxint1.physics.ox.ac.uk. This can be very large (>10GB)!'
        call('ssh vandenberg@cplxint1.physics.ox.ac.uk "tar -czpf - '+self.rootDIR+self.runfolder+'/*" | tar xzpf - -C ./',shell=True)
        print 'Download completed.\n'

    def tidyFulltree(self):
        print 'Cleaning up unneeded files ...'
        #call('cp -R .'+self.rootDIR+self.masterfolder+' '+self.masterfolder,shell=True)
        #call('rm -R ./home',shell=True)
        call('rm '+self.runfolder+'/scfly*',shell=True)
        call('rm '+self.runfolder+'/master*',shell=True)
        call('rm '+self.runfolder+'/i*/output/atomic*',shell=True)
        print 'sDone.\n'

    def cleanlocal(self):
        print 'Deleting all downloaded results from local tree...'
        for i in range (0,len(self.FEL_grid)):
            for j in range (0,self.INT_grid):

                folderpath = self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)+"/output"
                call('rm -R '+folderpath,shell=True)
        print 'Done.\n'

    def compact_spectra(self):
        print 'Binning and time-integrating spectra from local tree...'
        xaxis = np.linspace(self.wl_range[0],self.wl_range[-1],1000)

        for i in range (0,len(self.FEL_grid)):
            spectral_matrix = np.zeros([len(xaxis),self.INT_grid+1])
            spectral_matrix[0:,0] = xaxis;

            for j in range (0,self.INT_grid):
                folderpath = self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)+"/output"
                binned_spectrum = np.zeros(len(xaxis))

                # Unpack spectral tarball
                call('tar xjf '+folderpath+'/spec.tbz -C '+folderpath,shell=True)

                for k in range (2,self.TIM_grid+1):
                    spectrum = spec_files(folderpath,str(k).zfill(3)+'i'+str(i*int(self.INT_grid)+j+1))
                    binned_spectrum = np.add(binned_spectrum,spectrum.extract(xaxis))

                spectral_matrix[0:,j+1] = binned_spectrum
                # Clean up
                call('rm '+folderpath+'/sp*i*',shell=True)

            np.savetxt(self.runfolder+'/Time-integrated_sp'+str(int(self.FEL_grid[i]))+'eV.txt',spectral_matrix,delimiter='\t',fmt='%1.4e')
            print "Spectrum at "+str(self.FEL_grid[i])+" eV time-integrated."

        print 'Done.\n'

    def compact_spectra_local(self):
        print 'Binning and time-integrating spectra from local tree...'
        xaxis = np.linspace(self.wl_range[0],self.wl_range[-1],1000)

        for i in range (0,len(self.FEL_grid)):
            spectral_matrix = np.zeros([len(xaxis),self.INT_grid+1])
            spectral_matrix[0:,0] = xaxis;

            for j in range (0,self.INT_grid):
                folderpath = self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)
                binned_spectrum = np.zeros(len(xaxis))

                for k in range (2,self.TIM_grid+1):
                    spectrum = spec_files(folderpath,str(k).zfill(3)+'i'+str(i*int(self.INT_grid)+j+1))
                    binned_spectrum = np.add(binned_spectrum,spectrum.extract(xaxis))

                spectral_matrix[0:,j+1] = binned_spectrum

            np.savetxt(self.runfolder+'/Time-integrated_sp'+str(int(self.FEL_grid[i]))+'eV.txt',spectral_matrix,delimiter='\t',fmt='%1.4e')
            print "Spectrum at "+str(self.FEL_grid[i])+" eV time-integrated."

        print 'Done.\n'

    def compact_spectra_emissivity(self):
        print 'Binning and time-integrating spectra using emissivity from local tree...'
        xaxis = np.linspace(self.wl_range[0],self.wl_range[-1],1000)

        for i in range (0,len(self.FEL_grid)):
            spectral_matrix = np.zeros([len(xaxis),self.INT_grid+1])
            spectral_matrix[0:,0] = xaxis;

            for j in range (0,self.INT_grid):
                folderpath = self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)+"/output"
                binned_spectrum = np.zeros(len(xaxis))

                # Unpack spectral tarball
                call('tar xjf '+folderpath+'/spec.tbz -C '+folderpath,shell=True)

                for k in range (2,self.TIM_grid+1):
                    spectrum = spec_files(folderpath,str(k).zfill(3)+'i'+str(i*int(self.INT_grid)+j+1))
                    binned_spectrum = np.add(binned_spectrum,spectrum.extract_emissivity(xaxis))

                spectral_matrix[0:,j+1] = binned_spectrum
                # Clean up
                call('rm '+folderpath+'/sp*i*',shell=True)

            np.savetxt(self.runfolder+'/Time-integrated_emissivity_sp'+str(int(self.FEL_grid[i]))+'eV.txt',spectral_matrix,delimiter='\t',fmt='%1.4e')
            print "Spectrum at "+str(self.FEL_grid[i])+" eV time-integrated."

        print 'Done.\n'

    def compact_opacity(self):
        print 'Binning and time-integrating opacity from local tree...'
        xaxis = np.linspace(self.wl_range[0],self.wl_range[-1],1000)

        for i in range (0,len(self.FEL_grid)):
            opacity_matrix = np.zeros([len(xaxis),self.INT_grid+1])
            opacity_matrix[0:,0] = xaxis;

            for j in range (0,self.INT_grid):
                folderpath = self.runfolder+"/i"+str(i*int(self.INT_grid)+j+1)+"/output"
                binned_opacity = np.zeros(len(xaxis))
                for k in range (2,self.TIM_grid+1):
                    opacity = opacity_files(folderpath,str(k).zfill(3)+'i'+str(i*int(self.INT_grid)+j+1))
                    binned_opacity = np.add(binned_opacity,opacity.extract(xaxis))

                opacity_matrix[0:,j+1] = binned_opacity
                # Clean up
                call('rm '+folderpath+'/op*i*',shell=True)

            np.savetxt(self.runfolder+'/Time-integrated_op'+str(int(self.FEL_grid[i]))+'eV.txt',opacity_matrix,delimiter='\t',fmt='%1.4e')
            print "Opacity at "+str(self.FEL_grid[i])+" eV time-integrated."

        print 'Done.\n'




class spec_files:

    def __init__ (self,root,index):
        self.filename = root+'/sp'+index

    def extract(self,x_axis):
        data = np.loadtxt(self.filename,skiprows=1,usecols=(0,1))
        interpolation_spec = interpolate.interp1d(data[0:,0],data[0:,1],kind= 'linear')
        return interpolation_spec(x_axis)

    def extract_emissivity(self,x_axis):
        data = np.loadtxt(self.filename,skiprows=1,usecols=(0,3))
        interpolation_spec = interpolate.interp1d(data[0:,0],data[0:,1],kind= 'linear')
        return interpolation_spec(x_axis)

    def weights(self,variation,Peak_intensity,int_grid):
        if variation=='gaussian':
            peso = np.loadtxt("weights.dat",skiprows=1,usecols=(2,3))
        elif variation=='fscan':
            peso = np.loadtxt("weights.dat",skiprows=1,usecols=(2,4))
        else:
            print "Spot-distribution option not supported."

        # Normalize weights
        peso[0:,0] = np.multiply(peso[0:,0],Peak_intensity)
        peso[0:,1] = np.divide(peso[0:,1],np.sum(peso[0:,1]))

        f_i = interpolate.interp1d(peso[0:,0],peso[0:,1])
        interpolation_weights = extrap1d(f_i)
        return interpolation_weights(int_grid)

class opacity_files:

    def __init__ (self,root,index):
        self.filename = root+'/op'+index

    def extract(self,x_axis):
        data = np.loadtxt(self.filename,skiprows=1,usecols=(0,1))
        interpolation_opacity = interpolate.interp1d(data[0:,0],data[0:,1],kind= 'linear')
        return interpolation_opacity(x_axis)



def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike
