#!/usr/bin/env python3.6
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
            return strtobool(input().lower())
        except ValueError:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')

class initFolderStructure:
    'Creates folders to store processed data.'

    def __init__(self,filepath):
        print (os.path.realpath(filepath)+"/processed_data")
        print ((os.path.isdir(os.path.realpath(filepath)+"/processed_data")), os.path.dirname(os.path.realpath(filepath))+"/processed_data")

        if (os.path.isdir(os.path.join(os.path.realpath(filepath),"processed_data"))==False):
            print ("Newly created /processed_data directory")
            call(["mkdir",os.path.realpath(filepath)+"/processed_data"])
        elif (user_yes_no_query("Overwrite old /processed_data folder?")==True):
            print ("Overwrote old /processed_data directory with new empty directory")
            call(["rm","-r",os.path.realpath(filepath)+"/processed_data"])
            call(["mkdir",os.path.realpath(filepath)+"/processed_data"])
        else:
            print("Folder /processed_data already exists and is left untouched.")

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
              print ("\nInput parameters read from file: ", self.basepath+"/"+file)
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

              if self.material == "C":
                  self.Z = 6;
              elif self.material == "Mg":
                  self.Z = 12;
              elif self.material == "F":
                  self.Z = 9;
              elif self.material == "Al":
                  self.Z = 13;
              elif self.material == "Si":
                  self.Z = 14;
              else:
                  print ("Add material to input routine to determine self.Z")

    def displayCount(self):
        print ("Total amount of input parameters: %d" % loadInput.fileCount)
        return

class spectra:
    'Routines for time-integration of spectra, f-scan weighing and smoothening.'
    # If no time boundaries are set the default is to integrate across all time steps

    def __init__(self,inputParameters, i_start=None,i_end=None,t_start=None,t_end=None):
        # print "input.path %s" % (inputParameters.basepath)
        self.basepath = inputParameters.basepath
        # print 'basepath ', self.basepath
        # folder and log file to store results
        if not os.path.isdir(self.basepath+"/processed_data/spectra"):
            call(["mkdir",os.path.join(self.basepath,"processed_data/spectra")])
        self.logfile = open(os.path.join(self.basepath,"processed_data/spectra","intspec.log"),'w')

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
        if (os.path.isdir(self.basepath+"/processed_data/spectra")\
          and (user_yes_no_query("Time-integrate spectra and store in /processed_data/spectra folder?")==True)):
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
                file_specout = open(os.path.join(self.basepath,"processed_data/spectra","time-integrated_i"+str(i)),'w')

                file_specout.write("%s\t%s\n" %('E [eV]','Intensity'))
                for k in range(len(self.freq_grid)):
                    file_specout.write("%1.2f %1.4e\n" %(self.freq_grid[k],intensity[k]))
                    # print ("%1.2f %1.4e\n" %(self.freq_grid[k],intensity[k]))
                file_specout.close()
            print ("Time-integrated spectra written to folder ",self.basepath+"/processed_data/spectra/")
        else:
            print ("Time-integrated spectra already exist in folder %s." %(self.basepath+"/processed_data/spectra"))

        return
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
            datapath = os.path.join(self.basepath,"processed_data/spectra","time-integrated_i"+str(int(i+1)))
            data = np.loadtxt(datapath, skiprows=1)
            x = data[:,0] # energy, eV
            y = data[:,1] # intensity
            # Plot for inspection
            # if i==10:
            # plt.plot(x, y)
            # plt.show()
            total_spectrum += surface_weights[i]*np.interp(self.freq_grid,x,y)

        # Write weights to file
        file_weights = open(os.path.join(self.basepath,"processed_data","fscanweights"),'w')
        file_weights.write("%s\t%s\t%s\n" %('Index','Intensity','Weight'))
        for k in range(len(surface_weights)):
            file_weights.write("%d %1.4e %1.4e\n" %(k+1, I_fscan_unnormalised[k], surface_weights[k]))

        # Write spectrum to file
        file_specout = open(os.path.join(self.basepath,"processed_data/spectra","time-integrated_total"),'w')
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
        return

    def broaden(self, i_folders, sigma, line_shape):
        # Somehow this results in negative spectra. Fix first before applying.
        print ("Fix broaden/smoothen bug first before using further.")
        if (i_folders[-1]-i_folders[0]>1 and len(i_folders==2)):
            i_folders=np.linspace(i_folders[0],i_folders[1],i_folders[-1]-i_folders[0]+1, dtype = 'int')

        if (isinstance(i_folders,(list, tuple, np.ndarray))==True):
            for i in i_folders:
                # Save copy of original spectrum and update log file
                if not os.path.isdir(self.basepath+"/processed_data/spectra/oldcopies"):
                    print ('Created directory for old copies of spectra')
                    call(["mkdir",self.basepath+"/processed_data/spectra/oldcopies"])
                call("cp "+self.basepath+"/processed_data/spectra/time-integrated_i"+str(int(i)) +' '+self.basepath+"/processed_data/spectra/oldcopies/time-integrated_i"+str(int(i))+'_'+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S").replace(' ','_'),shell=True)
                self.logfile.write('Broadened i'+str(int(i))+' spectrum by a %s line shape with width: %1.2f \n' %(line_shape, sigma))
                self.logfile.write('Copy of i'+str(int(i))+' spectrum made with timestamp '+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")+ "\n" %())

                # Read spectra from file
                datapath = os.path.join(self.basepath,"processed_data/spectra","time-integrated_i"+str(int(i)))
                data = np.loadtxt(datapath, skiprows=1)
                x = data[:,0] # energy, eV
                y = data[:,1] # intensity

                # Construct line shape
                L = max(x)-min(x)
                x_convolve = np.linspace(-L/2,L/2,len(x))

                if line_shape == 'GAUSS':
                    y_convolve = lambda x_convolve: np.exp(-x_convolve**2/(2*sigma**2))
                    norm = integrate.quad(y_convolve,-L/2,L/2)[0]
                else:
                    print ("Currently only supports Gaussian lineshape. Supply argument \'GAUSS\'")
                    return
                # Convolve original spectrum with line shape for smoothening
                y_convolve = np.exp(-x_convolve**2/(2*sigma**2))/norm
                y_smooth = np.convolve(y, y_convolve, mode='same')
                y_smooth = y_smooth*np.trapz(x,y)/np.trapz(x,y_smooth)
                # Write result to file
                file_specout = open(datapath,'w')
                file_specout.write("%s\t%s\n" %('E [eV]','Intensity'))
                for k in range(len(x)):
                    file_specout.write("%1.2f %1.4e\n" %(x[k],y_smooth[k]))
                file_specout.close()
                return
        elif (i_folders == 'TOTAL'):
            # Save copy of original spectrum and update log file
            if not os.path.exists(self.basepath+"/processed_data/spectra/oldcopies"):
                print ('Created directory for old copies of spectra')
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
                y_convolve = lambda x_convolve: np.exp(-x_convolve**2/(2*sigma**2))
                norm = integrate.quad(y_convolve,-L/2,L/2)[0]
            else:
                print ("Currently only supports Gaussian lineshape. Supply argument \'GAUSS\'")
            y_convolve = np.exp(-x_convolve**2/(2*sigma**2))/norm
            y_smooth = np.convolve(y, y_convolve, mode='same')
            y_smooth = y_smooth*np.trapz(x,y)/np.trapz(x,y_smooth)
            # Write result to file
            file_specout = open(datapath,'w')
            file_specout.write("%s\t%s\n" %('E [eV]','Intensity'))
            for k in range(len(x)):
                file_specout.write("%1.2f %1.4e\n" %(x[k],y_smooth[k]))
            file_specout.close()
            return
        else:
            print ("Either operate on an numpy vector of input directories, or the total f-scan weighted spectrum. Check input arguments.")
            return

class extract:
    'Extract populations, temperature, density, rates'

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
            self.t_start = 1
        if (t_end == None):
            self.t_end = inputParameters.size_timegrid

    def superconfiguration(self,Z,charge_range,state):
        # Return superconfiguration object {"configs": [....], "indices":[...]} from atomic input for use in rate and population extraction
        configobject = {}
        superconfigs = []
        indices = []
        # Compatibility checks
        if (max(charge_range)==Z):
            print ("Maximum charge range must be smaller than Z.")
            return
        if (charge_range[-1]-charge_range[0]>1 and len(charge_range)==2):
            charge_range = np.linspace((charge_range[0]),(charge_range[1])+1,(charge_range[1])-(charge_range[0])+2, dtype= 'int')
        else:
            charge_range = np.linspace(charge_range[0],charge_range[-1]+1,len(charge_range)+1,dtype=int)

        # Determine configuration
        atomic_data_file = os.path.join('/usr/local/lib/atom/','atomic.inp.'+str(Z).zfill(2))
        if (os.path.isfile(atomic_data_file)==False):
            print ("Atomic data file (atomic.inp.xx) not found. Add atomic data file to folder \\usr\\local\\lib\\atom.")
            return
        periodic_table = {"1":"h_", "2":"he","3":"li","4":"be","5":"b_","6":"c_","7":"n_","8":"o_","9":"f_",\
            "10":"ne","11": "na", "12":"mg","13":"al","14":"si","15":"p_","16":"s_","17":"cl","18":"ar","19":"k_",\
            "20":"ca","21": "sc", "22":"ti","23":"v_","24":"cr","25":"mn","26":"fe","27":"co","28":"ni","29":"cu",\
            "30":"zn","31": "ga", "32":"ge","33":"as","34":"se","35":"br","36":"kr","37":"rb","38":"sr","39":"y_"}
        try:
            periodic_table[str(Z)]
        except KeyError:
            print ("Element isn't added to table yet. Update periodic table stored in \'superconfiguration\' function.")
            return
        num_electrons = np.full(charge_range.shape, Z, dtype='int') - charge_range #array of num electrons
        states = ['gs','single_ch','double_ch'] # non-core hole, single core-hole, double core-hole
        if state == states[0]:
            K_shell = 2
        elif state == states[1]:
            K_shell = 1
        elif state == states[2]:
            K_shell = 0
        else:
            print ("State is not supported. Choose from \'gs\',\'single_ch\',\'double_ch\'.")
        if (Z-max(num_electrons) - K_shell > 8):
            print ("Material not supported yet. Configurations can only be made for num_electrons =< 10.")
            return
        # Build superconfiguration strings
        for i in range(len(charge_range)):
            superconfigs.append(periodic_table[str(num_electrons[i])]+str(K_shell)+str(num_electrons[i]-K_shell))
        # Match superconfigs with atomic data strings to get full configuration and indices
        with open(atomic_data_file,'r') as atomfile:
            for line in atomfile:
                for i, config in enumerate(superconfigs):
                    if re.search(config,line):
                        superconfigs[i] = ' '+line.split()[3]
                        indices.append(line.split()[0].rjust(2)+'   '+line.split()[1].rjust(2))
        del superconfigs[-1]
        indices.reverse() # descend in charge state

        configobject["configs"] = superconfigs
        configobject["indices"] = indices
        return configobject

    def temperature_density(self,inputParameterse):
        # Set up directory
        if ((not os.path.isdir(self.basepath+"/processed_data/conditions")) or (user_yes_no_query("Overwrite /processed_data/conditions extracted plasma conditions folder?") == True)):
            call(["rm","-r",os.path.join(self.basepath,"processed_data/conditions")])
            call(["mkdir",self.basepath+"/processed_data/conditions"])
            print ("New folder created: %s" %(os.path.join(self.basepath,"processed_data/conditions")))
        else:
            print ("Folder already existed to store conditions: %s" %(os.path.join(self.basepath,"processed_data/conditions")))
            return
        self.conditions_logfile = open(self.basepath+"/processed_data/conditions/conditions.log",'w')
        self.conditions_logfile.write('Intensity folders considered for temperature-density evolution:\ti%s until i%s\n' % (str(self.i_start),str(self.i_end)))
        self.conditions_logfile.write('Time steps considered:\t%s until %s\n' %(str(self.t_start),str(self.t_end)))

        # Weigh conditions to f-scan weighted value
        if (os.path.isfile(self.basepath+"/processed_data/fscanweights")==False):
            print ("Temperature-density conditions can't be f-scan weighted as weights don't exist. Run spec.fscan first.")
        fscan_flag = ((os.path.isfile(self.basepath+"/processed_data/fscanweights")==True) and ((os.path.isfile(self.basepath+"/processed_data/conditions/trho_fscanweighted.txt")==False) or (user_yes_no_query("Overwrite f-scan weighted T-rho?")==True)))
        if (fscan_flag==True):
            data = np.loadtxt(self.basepath+"/processed_data/fscanweights", skiprows=1)
            weights = data[:,2]
            trhofscan_out = open(self.basepath+"/processed_data/conditions/trho_fscanweighted.txt",'w')
            trhofscan_out.write("%s\t%s\t%s\t%s\n" %('Time index','Time [s]','Temperature [eV]','Density [/cc]'))

        for i in range(self.i_start,self.i_end+1):
            i_folderpath = os.path.join(self.basepath, "i" + str(i),"output","zb.i"+str(i))
            data = np.loadtxt(i_folderpath, skiprows=1)
            trho_out = open(self.basepath+"/processed_data/conditions/trho_i%d.txt" %(i),'w')
            trho_out.write("%s\t%s\t%s\t%s\n" %('Time index','Time [s]','Temperature [eV]','Density [/cc]'))
            for j in range(data.shape[0]):
                trho_out.write("%d\t%1.4e\t%1.4f\t%1.4e\n" %(j+1,data[j,1],data[j,2],data[j,3])) # index, time, temperature, density
            trho_out.close()
            if fscan_flag==True:
                if (i==self.i_start):
                    T_fscan = np.zeros(data.shape[0])
                    rho_fscan = np.zeros(data.shape[0])
                # Add contribution to total f-scan
                T_fscan += weights[i-1]*data[:,2]/sum(weights[self.i_start:self.i_end])
                rho_fscan += weights[i-1]*data[:,3]/sum(weights[self.i_start:self.i_end])
                # Write out data per i folder

                if (i==self.i_end):
                    for j in range(data.shape[0]):
                        trhofscan_out.write("%d\t%1.4e\t%1.4f\t%1.4e\n" %(int(j+1),data[j,1],T_fscan[j],rho_fscan[j]))
                    trhofscan_out.close()
        return

    def populations(self, inputParameters, charge_range, state): # self, input_parameters, [charge start, charge end], state (gs, sch, dch)
        if (charge_range[-1]-charge_range[0]>1 and len(charge_range)==2):
            charge_range = np.linspace((charge_range[0]),(charge_range[1]),(charge_range[1])-(charge_range[0])+1, dtype= 'int')
        # Set up directory
        if ((not os.path.isdir(self.basepath+"/processed_data/populations")) or (user_yes_no_query("Overwrite /processed_data/populations extracted populations folder?") == True)):
            call(["rm","-r",os.path.join(self.basepath,"processed_data/populations")])
            call(["mkdir",self.basepath+"/processed_data/populations"])
            print ("New folder created: %s" %(os.path.join(self.basepath,"processed_data/populations")))
            print ('Populations extracted on '+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")+ "\n")
        else:
            print ("Folder already existed to store populations: %s" %(os.path.join(self.basepath,"processed_data/populations")))
            return
            exit()
        self.populations_logfile = open(self.basepath+"/processed_data/populations/populations.log",'w')
        self.populations_logfile.write('Intensity folders considered for population evolution:\ti%s until i%s\n' % (str(self.i_start),str(self.i_end)))
        self.populations_logfile.write('Time steps considered:\t%s until %s\n' %(str(self.t_start),str(self.t_end)))

        # Initialise parameters for f-scan weighted populations
        if (os.path.isfile(os.path.join(self.basepath,"processed_data/fscanweights"))==False):
            print ("Populations can't be f-scan weighted as weights don't exist. Run spec.fscan first.")
        fscan_flag = ((os.path.isfile(self.basepath+"/processed_data/fscanweights")==True) and\
            ((not os.path.isdir(os.path.join(self.basepath,"processed_data/populations", "z"+str(inputParameters.Z)+"_"+"_fscan")))\
             or (user_yes_no_query("Write out f-scan weighted populations?")==True)))
        if (fscan_flag==True):
            data = np.loadtxt(os.path.join(self.basepath,"processed_data/fscanweights"), skiprows=1)
            weights = data[:,2]
            pop_fscan = np.zeros([int(self.t_end)-int(self.t_start)+1, len(charge_range)]) # time, charge, per i-folder

            popfscan_file = open(os.path.join(self.basepath,"processed_data/populations", "z"+str(inputParameters.Z)+"_"+"_fscan"),'wb')
            self.populations_logfile.write('\nCalculated averaged populations from f-scan weighing.\n')
            self.populations_logfile.write("f-scan weights used:\n")
            self.populations_logfile.write('\t'.join(map(str,np.around(weights,5))))
            self.populations_logfile.write("\n")

        # Construct charge and material dependence to build population string in dedicated function
        popstring = self.superconfiguration(inputParameters.Z,charge_range,state)["configs"]
        print (popstring)
        # popstring = [" f_180002", " o_170002", " n_160002", " c_150002", " b_140002", " be130002", " li120002", " he110002"]

        # Loop over states
        states = ['gs','single_ch','double_ch'] # non-core hole, single core-hole, double core-hole
        for i in range(self.i_start,self.i_end+1):
            pop_out = np.zeros([int(self.t_end)-int(self.t_start)+1, len(popstring)]) # time, charge, per i-folder
            time_out = np.linspace(0.0,inputParameters.Tmax*1e-15,int(self.t_end)-int(self.t_start)+1)
            inputfile = os.path.join(self.basepath,"i"+str(i),"output",str(inputParameters.Z).zfill(2)+".i"+str(i))
            pop_outdirectory = os.path.join(self.basepath,"processed_data/populations")

            for k in range(len(popstring)):
                popdata = [];
                with open(inputfile) as popFile:
                    for line in popFile:
                        if re.match(popstring[k], line):
                            popline = line.split()[1:11]
                            popdata = popdata + popline
                        # elif re.match(' Time', line):
                        #     timeline = line.split()[1:11]
                        #     timedata = timedata + timeline
                pop_out[:,k]=popdata

            # Write out populations
            output_file = os.path.join(pop_outdirectory,"z"+str(inputParameters.Z)+"_"+state+"_i"+str(i))
            np.savetxt(output_file, np.concatenate((time_out[:,None],pop_out),axis=1), fmt='%1.4e', delimiter='\t', \
                newline='\n', header="%s\t\t%s\n%s\t\t%d%s%d\n%s\t\t%s\n%s\t%s" %("Material:",inputParameters.material,\
                "Charge range:",int(min(charge_range))," till ", int(max(charge_range)),"State:",state,"Time","\t".join(popstring)))

            if (fscan_flag==True):
                pop_fscan += weights[i-1]*pop_out
                if (i==self.i_end):
                    np.savetxt(popfscan_file, np.concatenate((time_out[:,None],pop_fscan),axis=1), fmt='%1.4e', delimiter='\t', newline='\n',\
                        header="%s\t\t%s\n%s\t\t%d%s%d\n%s\t\t%s\n%s\t%s" %("Material:",inputParameters.material,"Charge range:",int(min(charge_range)),\
                        " till ", int(max(charge_range)),"State:",state,"Time","\t".join(popstring)))
        return

    def rates(self,inputParameters,charge_range,state,rate_process): # self, input_parameters, [charge start, charge end], states element, rate_processes element
        np.set_printoptions(precision=4)
        ## Create directory and update log file
        if (charge_range[-1]-charge_range[0]>1 and len(charge_range)==2):
            charge_range = np.linspace((charge_range[0]),(charge_range[1]),(charge_range[1])-(charge_range[0])+1, dtype= 'int')
        rates_outdirectory = os.path.join(self.basepath,"processed_data/rates")

        if ((not os.path.isdir(rates_outdirectory)) or (user_yes_no_query("Overwrite /processed_data/rates extracted rates folder?") == True)):
            call(["rm","-r",rates_outdirectory])
            call(["mkdir",rates_outdirectory])
            print ("New folder created: %s" %(rates_outdirectory))
            print ('Rates extracted on '+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")+ "\n")
        else:
            print ("Folder already existed to store rates: %s" %(rates_outdirectory))
            return
        self.rates_logfile = open(os.path.join(rates_outdirectory,"rates.log"),'w')
        self.rates_logfile.write('Rates extracted on '+datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S")+ "\n")
        self.rates_logfile.write("Drawn from directory %s\n" %(self.basepath))
        self.rates_logfile.write('Intensity folders considered for extracted rates evolution:\ti%s until i%s\n' % (str(self.i_start),str(self.i_end)))
        self.rates_logfile.write('Time steps considered:\t%s until %s\n' %(str(self.t_start),str(self.t_end)))

        # Initialise parameters for f-scan weighted rates
        if (os.path.isfile(self.basepath+"/processed_data/fscanweights")==False):
            print ("Rates can't be f-scan weighted as weights don't exist. Run spec.fscan first.")
        fscan_flag = ((os.path.isfile(self.basepath+"/processed_data/fscanweights")==True) and\
            ((not os.path.isdir(os.path.join(self.basepath,"processed_data/rates", "z"+str(inputParameters.Z)+"_"+rate_process+"_fscan")))\
             or (user_yes_no_query("Write out f-scan weighted populations?")==True)))
        if (fscan_flag==True):
            data = np.loadtxt(self.basepath+"/processed_data/fscanweights", skiprows=1)
            weights = data[:,2]
            rates_fscan = np.zeros([int(self.t_end)-int(self.t_start), len(charge_range)]) # time, charge, per i-folder

            ratesfscan_file = open(os.path.join(self.basepath,"processed_data/rates", "z"+str(inputParameters.Z)+"_"+rate_process+"_fscan"),'wb')
            self.rates_logfile.write('\nCalculated averaged rates from f-scan weighted rates.\n')
            self.rates_logfile.write("f-scan weights used:\n")
            self.rates_logfile.write('\t'.join(map(str,np.around(weights,5))))
            self.rates_logfile.write("\n")

        # Loop over states and processes
        rate_processes = ['coll_ion','3body','auger'] # Collisional ionisation, 3-body recombination, auger decay
        states = ['gs','single_ch','double_ch'] # non-core hole, single core-hole, double core-hole
        ## Replace popstring and rate_idx by inherited function in 'extract' class
        popstring = self.superconfiguration(inputParameters.Z,charge_range,state)["configs"]
        rate_idx = self.superconfiguration(inputParameters.Z,charge_range,state)["indices"]
        # popstring = [" f_180002", " o_170002", " n_160002", " c_150002", " b_140002", " be130002", " li120002", " he110002"]
        # rate_idx = [" 9   12"," 8   12"," 7   12"," 6   10"," 5    8"," 4    6"," 3    4"," 2    2"," 1    1"]
        rate_idx_infile = 0
        if rate_process is rate_processes[0]:
            rate_idx_infile = 10
        elif  rate_process is rate_processes[1]:
            rate_idx_infile = 11
        elif rate_idx_infile is rate_processes[2]:
            rate_idx_infile = 9
            print("Warning: Make sure indices for Auger decay are properly implemented.")
        else:
            print ("Rate process is not supported yet. Choose between \'coll_ion\', \'3body\' and \'auger\'.")

        for i in range(self.i_start,self.i_end+1):
            rates_out = np.zeros([int(self.t_end)-int(self.t_start), len(popstring)]) # time, charge, per i-folder
            time_out = np.zeros(int(self.t_end)-int(self.t_start))
            # Extract rates from scfly output
            for j in range(int(self.t_start)+1,int(self.t_end)+1):
                rates_in = os.path.join(self.basepath,"i"+str(i),"output")+"/rate.%03d" %(j)
                # print (rates_in, rate_idx_infile) # correct file is being opened
                with open(rates_in) as rateFile:
                    time_out[j-1-int(self.t_start)] = rateFile.readline().split()[0]
                    for line in rateFile:
                        for k, pop in enumerate(popstring):
                            if re.match("i "+str(inputParameters.Z).zfill(2)+" "+rate_idx[k]+" "+rate_idx[k+1], line):
                                rates_out[j-1-int(self.t_start),k] = line.split()[rate_idx_infile]
                                # print("Extracted rate: ",line, line.split()[rate_idx_infile])
            # Write rates to file
            output_file = os.path.join(rates_outdirectory,"z"+str(inputParameters.Z)+"_"+rate_process+"_i"+str(i))
            np.savetxt(output_file, np.concatenate((time_out[:,None],rates_out),axis=1), fmt='%1.4e', delimiter='\t', newline='\n',\
                header="%s\t\t%s\n%s\t\t%d%s%d\n%s\t\t%s\n%s\t\t%s\n%s\t%s" %("Material:",inputParameters.material,"Charge range:",int(min(charge_range)),\
                " till ", int(max(charge_range)),"Rate process:",rate_process,"State:",state,"Time","\t".join(popstring)))
            # Compute f-scan weighted rates
            if (fscan_flag==True):
                rates_fscan += weights[i-1]*rates_out
                if (i == self.i_end):
                    np.savetxt(ratesfscan_file, np.concatenate((time_out[:,None],rates_fscan),axis=1), fmt='%1.4e', delimiter='\t', newline='\n',\
                        header="%s\t\t%s\n%s\t\t%d%s%d\n%s\t\t%s\n%s\t\t%s\n%s\t%s" %("Material:",inputParameters.material,"Charge range:",int(min(charge_range)),\
                        " till ", int(max(charge_range)),"Rate process:",rate_process,"State:",state,"Time","\t".join(popstring)))
        return





















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
                print ("Warning: Incorrect charge state.\nPlease choose a charge state between 0 and 5 for Carbon!")


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
                print ("Warning: Incorrect charge state.\nPlease choose a charge state between 0 and 11 for Mg!")


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
                print ("Warning: Incorrect charge state.\nPlease choose a charge state between 0 and 12 for Al!")


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
                print ("Warning: Incorrect charge state.\nPlease choose a charge state between 4 and 25 for Fe!")


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
                print ("Warning: Incorrect charge state.\nPlease choose a charge state between 2 and 28 for Cu!")
