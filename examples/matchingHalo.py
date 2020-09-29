# matchingHalo.py

import eagle
import numpy as N
import math
from scipy import stats
from time import time
import sys
 
# Simulation and snapshot 
sim_EA = '/cosma5/data/Eagle/ScienceRuns/Planck1/L0050N0752/PE/Z0p10_W1p00_E_3p0_0p3_ALPHA1p0e6_rhogas1_reposlim3p0soft_100mgas/data/'
sim_DM = '/cosma5/data/Eagle/ScienceRuns/Planck1/L0050N0752/PE/DMONLY_good_comp/data'
tag = '028_z000p000'
 
# Limits
fracToFind = 0.5
IDsToMatch = 50
massLimit = 2e9
N_min = 5000
 
mass_diff = 6       # mass ratio
max_distance = 8000 # kpc
 
 
# Constants
G = 4.302e4; # kpc (1e10 Mo)^-1 km^2 s^-2
unitMass = 1e10
unitLength = 1e3
Mpc = 3.08567758e22
kpc = 3.08567758e19
pc = 3.08567758e16
M_sun = 1.98855e30
 
numHaloes = 0
 
rank = int(sys.argv[1])
jobs = int(sys.argv[2])
 
for s in range(2):
 
    if s == 0:
        sim = sim_EA
    else:
        sim = sim_DM
 
    # General information
    numGroups = eagle.readAttribute("FOF", sim, tag, "/Header/TotNgroups")
    numSubGroups = eagle.readAttribute("SUBFIND", sim, tag, "/Header/TotNsubgroups")
    boxSize = eagle.readAttribute("PARTDATA", sim, tag, "/Header/BoxSize")
    hubbleParam = eagle.readAttribute("PARTDATA", sim, tag, "/Header/HubbleParam")
    H = eagle.readAttribute("FOF", sim, tag, "/Header/H(z)") * Mpc / 1000
    rho_crit = 3.*(H / unitLength)**2  / (8. * math.pi * G) * unitMass
    rho_bar = eagle.readAttribute("PARTDATA", sim, tag, "/Header/Omega0") * rho_crit
    redshift = eagle.readAttribute("PARTDATA", sim, tag, "/Header/Redshift")
    z_int = math.floor(redshift)
    z_dec = math.floor(10.*(redshift - z_int))
    expansionFactor = eagle.readAttribute("PARTDATA", sim, tag, "/Header/ExpansionFactor")
    physicalBoxSize = boxSize / hubbleParam
    softeningComoving = eagle.readAttribute("PARTDATA", sim, tag, "/RuntimePars/SofteningGas") * unitLength
    softeningMaxPhys = eagle.readAttribute("PARTDATA", sim, tag, "/RuntimePars/SofteningGasMaxPhys") * unitLength
    softening = min(softeningMaxPhys, softeningComoving * expansionFactor)
 
    # Print general info
    print("Sim:", sim))
    print("Redshift:" ,  redshift)
    print("Expansion factor:" ,  expansionFactor)
    print("H(z) =" ,  H )
    print("H0 =" ,  hubbleParam*100 )
    print("Box Size =" , physicalBoxSize , "Mpc.")
    print("Found" , numGroups , "haloes")
    print("Found" , numSubGroups , "sub-haloes")
    print("rho_crit =" , rho_crit , "M_sun / kpc^3")
    print("rho_bar =" , rho_bar , "M_sun / kpc^3"   )
    print("softening =" , softening , "kpc (", softeningComoving * expansionFactor, "," , softeningMaxPhys , ")")
    print("----------------------------------------------------------------")
 
    if numSubGroups > numHaloes:
        numHaloes = numSubGroups
 
 
particleIDs_EA = eagle.readArray("SUBFIND_PARTICLES", sim_EA, tag, "IDs/ParticleID", numThreads=16, oldSubfindFix=True)
particleIDs_DM = eagle.readArray("SUBFIND_PARTICLES", sim_DM, tag, "IDs/ParticleID", numThreads=16, oldSubfindFix=True)
 
 
M_EA = eagle.readArray("SUBFIND", sim_EA, tag, "Subhalo/Mass") * unitMass
CoP_EA = eagle.readArray("SUBFIND", sim_EA, tag, "Subhalo/CentreOfPotential") * unitLength
Grp_EA = eagle.readArray("SUBFIND", sim_EA, tag, "Subhalo/GroupNumber")
Sub_EA = eagle.readArray("SUBFIND", sim_EA, tag, "Subhalo/SubGroupNumber")
 
M_DM = eagle.readArray("SUBFIND", sim_DM, tag, "Subhalo/Mass") * unitMass
CoP_DM = eagle.readArray("SUBFIND", sim_DM, tag, "Subhalo/CentreOfPotential") * unitLength
Grp_DM = eagle.readArray("SUBFIND", sim_DM, tag, "Subhalo/GroupNumber")
Sub_DM = eagle.readArray("SUBFIND", sim_DM, tag, "Subhalo/SubGroupNumber")
 
length_EA = eagle.readArray("SUBFIND", sim_EA, tag, "Subhalo/SubLength")
lengthType_EA = eagle.readArray("SUBFIND", sim_EA, tag, "Subhalo/SubLengthType")
offset_EA = eagle.readArray("SUBFIND", sim_EA, tag, "Subhalo/SubOffset")
 
M_type_EA = eagle.readArray("SUBFIND", sim_EA, tag, "Subhalo/MassType") * unitMass
M_gas_EA = M_type_EA[:,0]
M_DM_EA = M_type_EA[:,1]
M_star_EA = M_type_EA[:,4]
 
length_DM = eagle.readArray("SUBFIND", sim_DM, tag, "Subhalo/SubLength")
lengthType_DM = eagle.readArray("SUBFIND", sim_DM, tag, "Subhalo/SubLengthType")
offset_DM = eagle.readArray("SUBFIND", sim_DM, tag, "Subhalo/SubOffset")
 
 
# Halo matched to the given ID
my_match_EA = N.zeros(numHaloes) - 1
my_match_DM = N.zeros(numHaloes) - 1
 
# ID that has been matched to the given halo
found_me_EA = N.zeros(numHaloes) - 1
found_me_DM = N.zeros(numHaloes) - 1
 
match_fraction_EA = N.zeros(numHaloes)
match_fraction_DM = N.zeros(numHaloes)
 
 
 
file = open("./matchedHalosSub_%d_z%dp%d_%d.dat"%(physicalBoxSize,z_int,z_dec, rank), 'w')
file.write("# All halos \n")
file.write("# box=%4.2fMpc z=%3.2f rho_crit=%5.3f epsilon=%4.3f IDsToMatch=%d\n"%(physicalBoxSize, redshift, rho_crit, softening, IDsToMatch))
file.write("# Gr_EAGLE    Sub_EAGLE    Gr_DMONLY    Sub_DMONLY    frac_1    frac_2    M_EAGLE [Msun]    M_DMONLY [Msun]    NDM_EAGLE     NDM_DMONLY\n")
 
 
# Helper function
def distance(x0, x1, dimensions):
    delta = N.abs(x0 - x1)
    delta = N.where(delta > 0.5 * dimensions, dimensions - delta, delta)
    return N.sqrt((delta ** 2).sum(axis=-1))
 
 
print("-----------------------------------------------------------")
print("Start loop over halos")
print("-----------------------------------------------------------")
 
# Loop over eagle halos
for i in range(rank, N.size(M_EA), jobs):
 
    sys.stdout.flush()
 
    # Consider only halos that are big enough
    if M_EA[i] > massLimit and lengthType_EA[i,1] > IDsToMatch:
        print("----\nFinding a match for halo (", Grp_EA[i], ",", Sub_EA[i], ") M=", M_EA[i], "CoP=", CoP_EA[i,:], "l_DM=", lengthType_EA[i,1])
 
        # Select 50 most bound particles of this halo
        halo_IDs = particleIDs_EA[long(offset_EA[i]): long(offset_EA[i])+long(length_EA[i])]
        mostBound_halo_IDs = halo_IDs[ halo_IDs % 2 == 0]
        mostBound_halo_IDs = mostBound_halo_IDs[0:IDsToMatch]
 
        # Loop over the other halo
        for j in range(N.size(M_DM)):
 
            # Consider only halos that are big enough
            if ( M_DM[j] > massLimit and lengthType_DM[j,1] > IDsToMatch ):
 
                # Eliminate stupid candidates
                if found_me_DM[j] != -1:
                    continue
 
                if( M_EA[i] > mass_diff * M_DM[j] or M_EA[i] < (1. / mass_diff) * M_DM[j] ):
                    continue
 
                if ( distance(CoP_EA[i], CoP_DM[j], N.ones(3)*physicalBoxSize * unitLength) > max_distance ):
                    continue
 
                # Select particles in this halo
                thisHalo_IDs = particleIDs_DM[long(offset_DM[j]): long(offset_DM[j]) + long(length_DM[j])]
 
                # Check whether the IDs from i are in j
                mask = N.in1d(thisHalo_IDs, mostBound_halo_IDs, assume_unique=True)
                matched_IDs = thisHalo_IDs[mask]
                count = N.size(matched_IDs)
 
                # Have we found enoguh particles ?
                if count >= fracToFind * IDsToMatch:
                    print("Matched halo (", Grp_EA[i], ",", Sub_EA[i], ") to halo (", Grp_DM[j], ",", Sub_DM[j], ") M=", M_DM[j], "CoP", CoP_DM[j,:], "l_DM=", lengthType_DM[j,1])
 
                    my_match_EA[i] = j;
                    match_fraction_EA[i] = (float)(count) / (IDsToMatch)
                    found_me_DM[j] = i;
 
                    print("Testing reversed match")
 
                    reversed_halo_IDs = particleIDs_DM[long(offset_DM[j]): long(offset_DM[j])+long(length_DM[j])]
                    reversed_halo_IDs = reversed_halo_IDs[0:IDsToMatch]
 
 
                    # Check whether the IDs from i are in j
                    reversed_mask = N.in1d(halo_IDs, reversed_halo_IDs, assume_unique=True)
                    reversed_matched_IDs = halo_IDs[reversed_mask]
                    reversed_count = N.size(reversed_matched_IDs)
 
                    if reversed_count >= fracToFind * IDsToMatch:
                        my_match_DM[j] = i;
                        match_fraction_DM[j] = (float)(reversed_count) / (IDsToMatch)
                        found_me_EA[i] = j;
                        print("Match confirmed. Fractions:", match_fraction_EA[i], match_fraction_DM[j])
 
 
                        # Print pair to file
                        file.write("%d %d %d %d %f %f %e %e %d %d\n"%( Grp_EA[i], Sub_EA[i], Grp_DM[j], Sub_DM[j], match_fraction_EA[i], match_fraction_DM[j], M_EA[i], M_DM[j],  lengthType_EA[i,1], lengthType_DM[j,1]))
                        file.flush()
 
 
                    else:
                        print("Match not confirmed.")
                    break
 
 
        if my_match_EA[i] == -1:
            print("No match found for halo", i)
 
    #else:
    #    print("----\nHalo too small (", Grp_EA[i], ",", Sub_EA[i], ") M=", M_EA[i], "CoP=", CoP_EA[i,:], "l_DM=", lengthType_EA[i,1])