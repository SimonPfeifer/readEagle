# haloProfile.py
# Plot the density, temperature and entropy profile of a halo
# Author: Matthieu Schaller
# Date: 11/02/2014

import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import eagle as E
from numpy import *


sim = "/cosma5/data/Eagle/ScienceRuns/Planck1/L0025N0376/PE/Z0p10_W1p00_E_3p0_0p3_ALPHA1p0e6_rhogas1_reposlim3p0soft_100mgas/data"
tag = "028_z000p000"

halo = 1 # Group number    

radialBinEdges = logspace(-4,0,30)
radialBins = 10**((log10(radialBinEdges[:-1]) + log10(radialBinEdges[1:])) / 2)

meanTemp = zeros(size(radialBinEdges)-1)
meanEntropy = zeros(size(radialBinEdges)-1)
meanDensity = zeros(size(radialBinEdges)-1)
medianTemp = zeros(size(radialBinEdges)-1)
medianEntropy = zeros(size(radialBinEdges)-1)
medianDensity = zeros(size(radialBinEdges)-1)


# Simulation information                                           
boxSize = E.readAttribute("SUBFIND", sim, tag, "/Header/BoxSize")
h = E.readAttribute("SUBFIND", sim, tag, "/Header/HubbleParam")
z = E.readAttribute("SUBFIND", sim, tag, "/Header/Redshift")

# Read in the data                                                  
pos = E.readArray("PARTDATA", sim, tag, "/PartType0/Coordinates")
rho = E.readArray("PARTDATA", sim, tag, "/PartType0/Density")
T = E.readArray("PARTDATA", sim, tag, "/PartType0/Temperature")
S = E.readArray("PARTDATA", sim, tag, "/PartType0/Entropy")
group = E.readArray("PARTDATA", sim, tag, "/PartType0/GroupNumber")

# Select particle from this halo                                   
pos = pos[group == halo, :]
rho = rho[group == halo]
T = T[group == halo]
S = S[group == halo]

# Get Halo center of potential                                     
CoP = E.readArray("SUBFIND_GROUP", sim, tag, "/FOF/GroupCentreOfPotential")[halo - 1,:]

# Compute radius                                                
r = sqrt((pos[:,0] - CoP[0])**2 + (pos[:,1] - CoP[1])**2 + (pos[:,2] - CoP[2])**2)  # Incorrect if CoP close to the boundary !!                                                    

# Sort the arrays by radius                                         
indices = argsort(r)
r = r[indices]
rho = rho[indices]
T = T[indices]
S = S[indices]

# Bin the data and compute mean                                     
meanDensity,b,c = stats.binned_statistic(r, rho, statistic='mean', bins=radialBinEdges)
meanTemp,b,c = stats.binned_statistic(r, T, statistic='mean', bins=radialBinEdges)
meanEntropy,b,c = stats.binned_statistic(r, S, statistic='mean', bins=radialBinEdges)

# Bin the data and compute median              
medianDensity,b,c = stats.binned_statistic(r, rho, statistic='median', bins=radialBinEdges)
medianTemp,b,c = stats.binned_statistic(r, T, statistic='median', bins=radialBinEdges)
medianEntropy,b,c = stats.binned_statistic(r, S, statistic='median', bins=radialBinEdges)

# Prepare entropy plot
params = {'axes.labelsize': 16, 'xtick.labelsize': 16, 'ytick.labelsize': 16, 'text.usetex': True, 'lines.linewidth' : 2}
rcParams.update(params)
figure(frameon=True)
subplot(111, xscale="log", yscale="log")

# Plot the simulation data                                
plot(radialBins, meanEntropy, 'r-', linewidth=2, label="Mean", markersize=6)
plot(radialBins, medianEntropy, 'b:', linewidth=2,label="Median",  markersize=6)

# Plot properties
legend(loc="upper left")
xlabel("r~[\\rm{kpc}]")
ylabel("S~[??]")
savefig("haloProfileEntropy.png")

# Prepare temperature plot
figure(frameon=True)
subplot(111, xscale="log", yscale="log")

# Plot the simulation data
plot(radialBins, meanTemp, 'r-', linewidth=2, label="Mean", markersize=6)
plot(radialBins, medianTemp, 'b:', linewidth=2, label="Median", markersize=6)

# Plot properties
legend(loc="upper left")
xlabel("r~[\\rm{kpc}]")
ylabel("T~[K]")
savefig("haloProfileTemp.png")

# Prepare density plot
figure(frameon=True)
subplot(111, xscale="log", yscale="log")

# Plot the simulation data
plot(radialBins, meanDensity, 'r-', linewidth=2, label="Mean", markersize=6)
plot(radialBins, medianDensity, 'b:',  linewidth=2, label="Median", markersize=6)

# Plot properties
legend(loc="lower left")
xlabel("r~[\\rm{kpc}]")
ylabel("$\\rho~[10^{10} M_\\odot \\rm{Mpc}^{-3}]$")

savefig("haloProfileDensity.png")