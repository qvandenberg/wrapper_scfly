## Plot spectra from SCFLY output
#  Written by Q.Y. van den Berg on 02.10.2017

import numpy as np
import math
from scipy import interpolate

## User input
root = '/Users/Quincy/Documents/SCFLY/8_ResCore_RT/12_IPD_experimental/longtime/1304eV_8e+16_';
expdata='/Users/Quincy/Documents/SCFLY/8_ResCore_RT/Expdata/Mg_highI/54nm/Mg_highI_all_data_calb_54nm.txt';
datasets = ['ekCF','spCF','ekLotz', 'ekBC'];
legendstring = ['Experimental','EK CF', 'SP CF', 'EK Lotz', 'EK BC'];


# SCFLY grid settings
Igrid = np.linspace(1,20,20);
time = [0, 1000, 300];
E = [1200, 1400, 500]; # Energy range in eV [E_min E_max #steps]
tsize = 5.40e-06; # Target thickness from hist file - 25 nm, currently not used

# f-scan parameters
fss = [0.957467, 3.54691, 0.46181, 0.042533, 0.1929, 0.21648]; # hyper-gaussian parameters
S = np.linspace(0,1000,5000); # spot surface grid

## no user input below this line

dt = (time[1]-time[0])/time[2];
Y_plot = np.zeros((len(datasets),E[2]))
Egrid = np.linspace(E[0],E[1],E[2]);
Totspec = []

for k in range(len(datasets)):
    for i in range(len(Igrid)):
        for j in range(2,time[2]):
            filename =[]
            if (j<10):
                filename =  root + datasets[k]+ '/i' + str(int(Igrid[i])) + '/output/out.00' + str(j)
            elif (j<100):
                filename =  root + datasets[k]+ '/i' + str(int(Igrid[i])) + '/output/out.0' + str(j)
            elif (j<1000):
                filename =  root + datasets[k]+ '/i' + str(int(Igrid[i])) + '/output/out.' + str(j)
            else:
                print "Amount of timesteps is too large. Try again for <1000 timesteps."
            # print "filename: ", filename
            if ((i==0) & (j==2)):
                 Totspec = np.zeros((len(Igrid),len(Egrid)));
                 print "Working on dataset: ", datasets[k]

            # Read data from filename
            tempspec = [];
            egridtemp = [];
            with open(filename) as specfile:
                for _ in xrange(2):
                    next(specfile)
                for line in specfile:
                    egridtemp.append(float(line.split()[0]));
                    tempspec.append(dt*float(line.split()[3]));

            # interpolate to match energy grid from plot
            f = interpolate.interp1d(egridtemp, tempspec, kind='linear',bounds_error=False,fill_value='extrapolate',assume_sorted="True") #, bounds_error="False", fill_value="None")
            # integrate in time and add result to time-integrated spectrum
            if (math.isnan(np.sum(f(Egrid)))==False):
                Totspec[i,:] += f(Egrid)
            Totspec[i,:] = map(lambda x: max(x,0), Totspec[i,:])

    # f-scan weighing
    Int_root = root + datasets[0]+'/intensity_grid.dat';
    Intgrid = [];
    with open(Int_root) as intfile:
        for line in intfile:
            Intgrid.append(float(line.split()[1]));
    # Both Intgrid and fS normalised to peak at 1
    Intgrid = map(lambda x: x / np.max(Intgrid), Intgrid)
    fS = fss[0]*np.exp(-(S/fss[1]))**fss[2] + fss[3]*np.exp(-(S/fss[4])) **fss[5];

    gw = np.zeros((len(Intgrid)));
    idx = np.zeros(len(Intgrid),dtype=np.int8);
    for i in range(0,len(Intgrid)):
        # Find index that maps Intgrid onto fS
        idx[i] = int((np.abs(fS-Intgrid[i])).argmin())
        if (i<(len(Intgrid)-1)):
            gw[i] = abs((S[idx[i+1]]-S[idx[i]]))
        elif (i ==(len(Intgrid)-1)):
            gw[i] = gw[i-1]
    # Normalise sum of f-scan weights to one
    gw = map(lambda x: x / np.sum(gw), gw)

    # Calculate total spectrum per dataset
    for i in range(0,len(Intgrid)):
        Totspec[i,:] *= gw[i]
        Y_plot[k,:] = np.add(Y_plot[k,:],Totspec[i,:])
        # Y_EK_direct(k,:) = Broaden(X(:),Y_EK_direct(k,:),1);

    # Write f-scan weighted and time integrated spectrum to output file
    file_out = open(root+datasets[k]+'/timespace_integrated_spectrum.out','w+a')
    file_out.write("Energy (eV)\tSpectrum\n")
    for i in range(len(Egrid)):
        file_out.write(str('{0:.2f}'.format(Egrid[i]))+"\t"+str('{0:.2f}'.format(Y_plot[k,i]))+"\n")



#
# %% Read in experimental spectrum
# expdata=importfile('/Users/Quincy/Documents/SCFLY/8_ResCore_RT/Expdata/Mg_highI/54nm/Mg_highI_all_data_calb_54nm.txt', 2, 2048);
# expdata = expdata(614:1000,1:2);
# expdata(:,2) = expdata(:,2)*1;
#
# %% Renormalize calculated spectra
#
# for i=1:length(dataset)%-1
# if i==1
# peaknorm = max(squeeze(expdata(:,2)));
# %peaknorm = 2800;
# end
# Y_EK_direct(i,:) = Y_EK_direct(i,:)*peaknorm/max(squeeze(Y_EK_direct(i,700:end)));
# if i==length(dataset)-1
# %    clear peaknorm
# end
# end
#
# %% Plot scaled spectra
# % Manually written to plot scaling {5, 3, 2.5, 1}
# % Remember to also (manually) renormalise central peak magnitude
# Eshift = .5; % shift energy axis
# figure,plot(squeeze(expdata(:,1))-1.5,expdata(:,2),'LineWidth',3)
# hold on
# %legendstring = ['Experimental',dataset];
# if length(dataset)>0
# plot(X+Eshift,squeeze(Y_EK_direct(1,:)),'-.','Linewidth',3) % Plots weighted spectrum. Index 2 = F_180, 9 = He110
# end
# if length(dataset)>1
# plot(X+Eshift,squeeze(Y_EK_direct(2,:)),'-.','Linewidth',3) % Plots weighted spectrum. Index 2 = F_180, 9 = He110
# end
# if length(dataset)>2
# plot(X+Eshift,squeeze(Y_EK_direct(3,:)),'-.','Linewidth',3) % Plots weighted spectrum. Index 2 = F_180, 9 = He110
# end
# if length(dataset)>3
# plot(X+Eshift,squeeze(Y_EK_direct(4,:)),'-.','Linewidth',3) % Plots weighted spectrum. Index 2 = F_180, 9 = He110
# end
# hold off
#
# xlabel('Energy, eV','Interpreter','latex','fontsize',25,'fontweight','bold')
# ylabel('Intensity a.u.' ,'Interpreter','latex','fontsize',25,'fontweight','bold')
# set(gca,'fontsize',22)
#
# legend(legendstring)
# h = legend;
# set(h,'Interpreter','latex','fontsize',25)''
# set(gca,'fontsize',25,'linewidth',2)
# set(gca, 'XTick', [1260 1280 1300 1320 1340 1360])
# set(gca, 'YTick', [0 1e3 2e3 3e3])
# axis([1250 1365 -1 3050])
#
# % save2pdf('/Users/Quincy/Documents/SCFLY/8_ResCore_RT/12_IPD_experimental/Figs/spIPD_4models.pdf',gcf,5000)
