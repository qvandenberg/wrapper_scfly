##  Written by Muhamad on 18.10.2017 to plot scfly spectra directly from directory
#
# Call with "python plot_spec.py i1/output" or wherever the output
# folder exact path is at.
#

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plotting the intensity given the output directory names")
parser.add_argument("outdir", type=str, help="Output directory to be plotted")
parser.add_argument("--intensity", dest="is_intensity", action="store_const", const=True, default=False, help="If specified, plot the intensity, otherwise, plot the emissivity")
parser.add_argument("--fortran", dest="is_fortran", action="store_const", const=True, default=False, help="If specified, the output file is from Fortran spectral code")

args = parser.parse_args()

# frequency points
freq_pts = np.linspace(1220, 1370, 2000)

# where the total intensity will be stored
intensity = np.zeros_like(freq_pts)

# iterate for every time step
for i in range(2, 51):

    # read the output
    if args.is_fortran:
        data = np.loadtxt(os.path.join(args.outdir, "sp%03d%s" % (i, "i1")), skiprows=1)
        freqs = data[:,0]
        idx = 1 if args.is_intensity else 3
        intens = data[:,idx]
    else:
        data = np.loadtxt(os.path.join(args.outdir, "out.%03d" % i), skiprows=2)
        freqs = data[:,0]
        idx = 3 if args.is_intensity else 1 # idx = 3 for intensity and 1 for emissivity
        intens = data[:,idx] 

    # interpolate the intensity
    intensity += np.interp(freq_pts, freqs, intens)

# Write out intensities to file
filename = "Extracted/time_int_spectrum"
file_specout = open(filename,'w+a')
for i in range(len(freq_pts)):
    file_specout.write("%1.2f %1.2f\n" %(freq_pts[i],intensity[i]))
    
#    file_specout.write(str(freq_pts[i])+"\t"+str(intensity[i])+"\n")
print "len freq len int ",len(freq_pts),len(intensity)
print "Written spectrum to file: "+filename
plt.plot(freq_pts, intensity)
plt.show()
