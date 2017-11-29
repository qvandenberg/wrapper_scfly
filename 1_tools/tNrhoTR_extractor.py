#!/usr/bin/env python2.7

#__________________________________________________________________________
#
#  Created by Sam M. Vinko on 15/09/2015.
#  QvdB: Added recombination & loop over charge states on 09/05/2016.
#__________________________________________________________________________

# operate with python rate_extractor.py i1 i20
import sys, re

# path to scfly output
root   = "/Users/Quincy/Documents/Code/SCFLY_analysis/20171018_SPipd_old/1304eV_8e+16_ekLotz"

# Time steps and atomic number
start_i = 2
finis_i = 51
z = "12"
# Species to extract population and rates from. Last name not included
popstring = [ " b_140002"]
rate_idx = [" 5    8"," 4    6"]
#popstring = [" f_180002", " o_170002", " n_160002", " c_150002", " b_140002", " be130002", " li120002", " he110002"]
#rate_idx = [" 9   12"," 8   12"," 7   12"," 6   10"," 5    8"," 4    6"," 3    4"," 2    2"," 1    1"]

## No user input under this line

folder_beg = sys.argv[1]
folder_end = sys.argv[2]

for k in range(len(popstring)):
#    filename_outmaster = root + "/Extracted/master_"+root[-4:]+str(k+1)+".dat";
    filename_outmaster = root + "/Extracted/master_"+root.split('_')[-1]+str(k+1)+".dat";

    
    filename_outP = root + "/Extracted/N"+str(k+1)+".dat";
    filename_outRi =root + "/Extracted/R"+str(k+1)+".dat";
#    filename_outP = "/Users/Quincy/Documents/SCFLY/8_ResCore_RT/Extracted/Lotz/N"+str(k+1)+"_"+folder_beg+".dat";
#    filename_outRi = "/Users/Quincy/Documents/SCFLY/8_ResCore_RT/Extracted/Lotz/Ri"+str(k+1)+"_"+folder_beg+".dat";
#    filename_outRr = "/Users/Quincy/Documents/SCFLY/8_ResCore_RT/Extracted/Lotz/Rr"+str(k+1)+"_"+folder_beg+".dat";
    file_outmaster = open(filename_outmaster,'w+a')
#    file_outP = open(filename_outP,'w+a')
#    file_outRr = open(filename_outRr,'w+a')
#    file_outRi = open(filename_outRi,'w+a')

    for intensity in range(int(folder_beg[1:]),int(folder_end[1:])+1):
        print "Working on i"+str(intensity)+" ..."
        
        # Read from population file
        popdata_150 = [];
        filename_P = root+"/i"+str(intensity)+"/12.i"+str(intensity)
        with open(filename_P) as popFile:
            for line in popFile:
              if re.match(popstring[k], line):
                   popline = line.split()[1:11]
                   popdata_150 = popdata_150 + popline
    
        # Read from rate file
        ion_rate =[];
        # Open rate file for ionisation
        for i in range(start_i,finis_i):
            filename_Ri = root+"/i"+str(intensity)+"/rate."+ str(i).zfill(3)
            with open(filename_Ri) as rateFile:
                for line in rateFile:
                    if re.match("i "+z+" "+rate_idx[k]+" "+rate_idx[k+1], line):
                      data = line.split()[10]
                      ion_rate.append(data)
        # Read from zb file for rho and T
        zbT = [];
        zbrho = [];
        zbt =[];
        filename_zb = root+"/i"+str(intensity)+"/zb.i"+str(intensity)
        for i in range(start_i,finis_i):
            with open(filename_zb) as zbFile:
                for line in zbFile:
                    if re.search(" "+str(i)+" ", line):
                        zbt.append(line.split()[1]) # time
                        zbT.append(line.split()[2]) # temperature
                        zbrho.append(line.split()[3]) # density

                                  
        # Write t, N, rho, T, Rate strings to output file
        file_outmaster.write("time(fs)\tN(cm-3)\trho_e\tT(eV)\tRate(cm-3s-1)\n")
        popdata_150.pop(0) # make equal length to other lists
        for i in range(1,len(popdata_150)):
            file_outmaster.write("%1.2e %1.2e %1.2e %1.2e %1.2e\n" %(float(zbt[i]),float(popdata_150[i]),float(zbrho[i]),float(zbT[i]),float(ion_rate[i])))

#        # Open rate file for recombination
#        for i in range(start_i,finis_i):
#            filename_Rr = root+"/i"+str(intensity)+"/rate."+ str(i).zfill(3)
#            with open(filename_Rr) as rateFile:
#                for line in rateFile:
#                    if re.match("i "+z+" "+rate_idx[k]+" "+rate_idx[k+1], line):
#                      data = line.split()[11]
#                      file_outRr.write(str(data)+"\t")
#        file_outRr.write("\n")

    print 'Written data to file: ', filename_outmaster

