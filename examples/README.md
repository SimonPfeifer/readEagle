# Examples

Below are compiled an (almost) exhaustive list of data types and outputs from a typical simulation along with a short description.

Recall that the general way to read in data using `readEagle` looks something like this:

```python
import eagle as E
 
attr = E.readAttribute(file_type,  directory,  tag,   attribute)
arr = E.readArray(file_type,  directory,  tag,   array)
```

## File type

The first argument required by `readEagle` is a `file_type`, a string describing the type of file and data to read. The allowed values of this argument are:  

|fileType | Description | Example of data that can be read |
| --- | --- | --- |
| FOF | FoF group informations | Group centre of mass, group length, group star formation rate |
| FOF_PARTICLES | IDs of the particles in a FOF group | Particle IDs |
| SNIP_FOF | FoF group informations (snipshot) | Group centre of mass, group length, group star formation rate |
| SNIP_FOF_PARTICLES | IDs of the particles in a FOF group (snipshot) | Particle IDs |
| PARTDATA | Particles that are in a FOF group | Particle Mass, velocity, entropy, stellar age |
| SNIP_PARTDATA | Particles that are in a FOF group (snipshot) | Particle Mass, velocity, entropy, stellar age |
| SNAPSHOT | Full information about all particles | Particle Mass, velocity, entropy, stellar age |
| SNIPSHOT | Reduced information about all particles | Mass, velocity |
| SUBFIND | Subhalo information | Subhalo mass, subhalo centre of potential |
| SUBFIND_GROUP | Subfind halo information | Group centre of potential, M_200, R_500 |
| SUBFIND_PARTICLES | IDs of the particles in a subhalo | Particle IDs |
| SNIP_SUBFIND | Subhalo information (snipshot) | Subhalo mass, subhalo centre of potential |
| SNIP_SUBFIND_GROUP | Subfind halo information (snipshot) | Group centre of potential, M_200, R_500 |
| SNIP_SUBFIND_PARTICLES | IDs of the particles in a subhalo (snipshot) | Particle IDs |

## Array data
The last argument required by `readEagle` refers directly to the HDF5 data within the given `fileType` and is written in a directory-like format:
```python
array_example = "/path/to/data"
```
You can use the lists below to find the exact data location. A helpful command is `h5ls` which you can use to browse the data structure of a HDF5 simulation output file. A quick example would be:
```
h5ls groups_032.0.hdf5/FOF
```
For a more user friendly way of browsing a HDF5 file you can use HDF5View which presents the file structure in a GUI.

### Snapshot
Sanpshots store the complete information about all particles in the simulation with associated particle type number:

| Number | Particle |
| --- | --- |
| 0 | Gas | 
| 1 | DM | 
| 2 | N/A | 
| 3 | N/A | 
| 4 | Stars |
| 5 | BHs |

The particle data available for each particle type is given in the table below.
| Particle Type | Array Name | Array Description |
| ---------- |----------------------------------------- | --- |
| PartType0  |     AExpMaximumTemperature               |Expansion factor a when particle had highest temperature|
| PartType0  |     Coordinates                          |Co-moving coordinates. Physical: r = a x = Coordinates h to the power of-1 a U_L [cm]|
| PartType0  |     Density                              |Co-moving mass densities. Physical rho = Densities h to the power of2 a to the power of-3 U_M U_L to the power of-3 [g/cm to the power of3]|
| PartType0  |     ElementAbundance/Carbon              | |
| PartType0  |     ElementAbundance/Helium              | |
| PartType0  |     ElementAbundance/Hydrogen            | |
| PartType0  |     ElementAbundance/Iron                | |
| PartType0  |     ElementAbundance/Magnesium           | |
| PartType0  |     ElementAbundance/Neon                | |
| PartType0  |     ElementAbundance/Nitrogen            | |
| PartType0  |     ElementAbundance/Oxygen              | |
| PartType0  |     ElementAbundance/Silicon             | |
| PartType0  |     Entropy                              |Particle entropy. Physical s = Entropy h to the power of(2-2*GAMMA) UnitPressure UnitDensity to the power of-GAMMA|
| PartType0  |     GroupNumber                          |FoF group number particle is in|
| PartType0  |     HostHalo_TVir_Mass                   |Estimate of halo's virial temperature, calculated from the DM halo mass.   T_vir = (MEANMOLIONIZED * PROTONMASS / 3. / BOLTZMANN) * (G * m200 * H(z)) to the power of(2./3.) [K]|
| PartType0  |     InternalEnergy                       |Thermal energy per unit mass. Physical u = InternalEnergy U_V to the power of2 [(cm/s) to the power of2]|
| PartType0  |     IronMassFracFromSNIa                 |Iron mass from SNIa divided by particle mass. (Initial particle mass for stars)|
| PartType0  |     Mass                                 |Particle mass. Physical m = Mass h to the power of-1 U_M [g]|
| PartType0  |     MaximumTemperature                   |Maximum temperature ever reached by a particle [K]|
| PartType0  |     MetalMassFracFromAGB                 |Metal mass from AGB and their progenitors divided by particle mass. (Initial particle mass for stars)|
| PartType0  |     MetalMassFracFromSNII                |Metal mass from SNII and their progenitors divided by particle mass. (Initial particle mass for stars)|
| PartType0  |     MetalMassFracFromSNIa                |Metal mass from SNIa divided by particle mass. (Initial particle mass for stars)|
| PartType0  |     MetalMassWeightedRedshift (A bug affects this quantity.)    |Metal mass weighted redshift at which particle was enriched.|
| PartType0  |     Metallicity                          |Mass fraction of elements heavier than Helium|
| PartType0  |     OnEquationOfState                    |Star-formation flag. 0 if has never been star-forming, +ve if currently sf, -ve if not currently sf, value indicates aexp at which it obtained its current state|
| PartType0  |     ParticleIDs                          |Unique particle identifier|
| PartType0  |     SmoothedElementAbundance/Carbon     | |
| PartType0  |     SmoothedElementAbundance/Helium     | |
| PartType0  |     SmoothedElementAbundance/Hydrogen   | |
| PartType0  |     SmoothedElementAbundance/Iron       | |
| PartType0  |     SmoothedElementAbundance/Magnesium  | |
| PartType0  |     SmoothedElementAbundance/Neon       | |
| PartType0  |     SmoothedElementAbundance/Nitrogen   | |
| PartType0  |     SmoothedElementAbundance/Oxygen     | |
| PartType0  |     SmoothedElementAbundance/Silicon    | |
| PartType0  |     SmoothedIronMassFracFromSNIa         |Smoothed mass from SNIa divided by particle mass. (Initial particle mass for stars)|
| PartType0  |     SmoothedMetallicity                  |Smoothed mass fraction of elements heavier than Helium|
| PartType0  |     SmoothingLength                      |Co-moving smoothing length. Physical h = SmoothingLength h to the power of-1 a U_L [cm]|
| PartType0  |     StarFormationRate                    |Gas star formation rate in solar masses / yr|
| PartType0  |     SubGroupNumber                       |Subgroup number particle is in|
| PartType0  |     Temperature                          |Temperature [K]|
| PartType0  |     TotalMassFromAGB                     |Total mass received from AGB and their progenitors|
| PartType0  |     TotalMassFromSNII                    |Total mass received from SNII and their progenitors|
| PartType0  |     TotalMassFromSNIa                    |Total mass received from SNIa|
| PartType0  |     Velocity                             |Co-moving velocities. Physical v_p = a dx/dt  = Velocities a to the power of1/2 U_V [cm/s]|
| PartType1  |     Coordinates                          |Co-moving coordinates. Physical: r = a x = Coordinates h to the power of-1 a U_L [cm]|
| PartType1  |     GroupNumber                          |FoF group number particle is in|
| PartType1  |     ParticleIDs                          |Unique particle identifier|
| PartType1  |     SubGroupNumber                       |Subgroup number particle is in|
| PartType1  |     Velocity                             |Co-moving velocities. Physical v_p = a dx/dt  = Velocities a to the power of1/2 U_V [cm/s]|
| PartType4  |      AExpMaximumTemperature              |Expansion factor a when particle had highest temperature|
| PartType4  |      BirthDensity                        |Local gas density (physical units)  when a star particle was born. No a-factor correction as the a-factor at birth time is factored in.|
| PartType4  |      Coordinates                         |Co-moving coordinates. Physical: r = a x = Coordinates h to the power of-1 a U_L [cm]|
| PartType4  |      ElementAbundance/Carbon            | |
| PartType4  |      ElementAbundance/Helium            | |
| PartType4  |      ElementAbundance/Hydrogen          | |
| PartType4  |      ElementAbundance/Iron              | |
| PartType4  |      ElementAbundance/Magnesium         | |
| PartType4  |      ElementAbundance/Neon              | |
| PartType4  |      ElementAbundance/Nitrogen          | |
| PartType4  |      ElementAbundance/Oxygen            | |
| PartType4  |      ElementAbundance/Silicon           | |
| PartType4  |      Feedback_EnergyFraction             |Energy fraction used for SNII feedback (no units).|
| PartType4  |      GroupNumber                         |FoF group number particle is in|
| PartType4  |      HostHalo_TVir                       |Halo's virial temperature used in Type II SNe feedback [K]|
| PartType4  |      HostHalo_TVir_Mass                  |Estimate of halo's virial temperature, calculated from the DM halo mass.   T_vir = (MEANMOLIONIZED * PROTONMASS / 3. / BOLTZMANN) * (G * m200 * H(z)) to the power of(2./3.) [K]|
| PartType4  |      InitialMass                         |Star particle mass at formation time. Physical m = InitialMass h to the power of-1 U_M [g]|
| PartType4  |      IronMassFracFromSNIa                |Iron mass from SNIa divided by particle mass. (Initial particle mass for stars)|
| PartType4  |      Mass                                |Particle mass. Physical m = Mass h to the power of-1 U_M [g]|
| PartType4  |      MaximumTemperature                  |Maximum temperature ever reached by a particle [K]|
| PartType4  |      MetalMassFracFromAGB                |Metal mass from AGB and their progenitors divided by particle mass. (Initial particle mass for stars)|
| PartType4  |      MetalMassFracFromSNII               |Metal mass from SNII and their progenitors divided by particle mass. (Initial particle mass for stars)|
| PartType4  |      MetalMassFracFromSNIa               |Metal mass from SNIa divided by particle mass. (Initial particle mass for stars)|
| PartType4  |      MetalMassWeightedRedshift (A bug affects this quantity.)           |Metal mass weighted redshift at which particle was enriched.|
| PartType4  |      Metallicity                         |Mass fraction of elements heavier than Helium|
| PartType4  |      ParticleIDs                         |Unique particle identifier|
| PartType4  |      PreviousStellarEnrichment           |This is the expansion factor when the star last did enrichment.|
| PartType4  |      SmoothedElementAbundance/Carbon    | |
| PartType4  |      SmoothedElementAbundance/Helium    | |
| PartType4  |      SmoothedElementAbundance/Hydrogen  | |
| PartType4  |      SmoothedElementAbundance/Iron      | |
| PartType4  |      SmoothedElementAbundance/Magnesium | |
| PartType4  |      SmoothedElementAbundance/Neon      | |
| PartType4  |      SmoothedElementAbundance/Nitrogen  | |
| PartType4  |      SmoothedElementAbundance/Oxygen    | |
| PartType4  |      SmoothedElementAbundance/Silicon   | |
| PartType4  |      SmoothedIronMassFracFromSNIa        |Smoothed mass from SNIa divided by particle mass. (Initial particle mass for stars)|
| PartType4  |      SmoothedMetallicity                 |Smoothed mass fraction of elements heavier than Helium|
| PartType4  |      SmoothingLength                     |Co-moving smoothing length. Physical h = SmoothingLength h to the power of-1 a U_L [cm]|
| PartType4  |      StellarEnrichmentCounter            |The counter shows the number of time steps since enrichment was last done.|
| PartType4  |      StellarFormationTime                |Expansion factor a when star particle was born|
| PartType4  |      SubGroupNumber                      |Subgroup number particle is in|
| PartType4  |      TotalMassFromAGB                    |Total mass received from AGB and their progenitors|
| PartType4  |      TotalMassFromSNII                   |Total mass received from SNII and their progenitors|
| PartType4  |      TotalMassFromSNIa                   |Total mass received from SNIa|
| PartType4  |      Velocity                            |Co-moving velocities. Physical v_p = a dx/dt  = Velocities a to the power of1/2 U_V [cm/s]|
| PartType5  |      BH_AccretionLength                  |BH smoothing length.|
| PartType5  |      BH_CumlAccrMass                     |Cumulative mass accreted by largest progenitor of this BH. Physical m = Mass h to the power of-1 U_M [g]|
| PartType5  |      BH_CumlNumSeeds                     |Cumulative number of BH seeds swallowed by this BH.|
| PartType5  |      BH_Density                          |Co-moving black hole densities. Physical rho = Densities h to the power of2 a to the power of-3 U_M U_L to the power of-3 [g/cm to the power of3]|
| PartType5  |      BH_EnergyReservoir                  |Black hole energy reservoir for thermal feedback.|
| PartType5  |      BH_FormationTime                    |Expansion factor a when BH particle was born|
| PartType5  |      BH_Mass                             |BH mass. Physical m = Mass h to the power of-1 U_M [g]|
| PartType5  |      BH_Mdot                             |BH accretion rate. Physical mdot = BH_Mdot h to the power of-1 U_M /U_T [g/s]|
| PartType5  |      BH_MostMassiveProgenitorID          |Unique ID of the most massive progenitor of this BH. At each merger event, the ID of the most massive of the two merging BHs is stored in this array.|
| PartType5  |      BH_Pressure                         |Black hole surrounding gas pressure. Physical P = Pressure h to the power of2 a to the power of(-3*GAMMA) U_M U_V to the power of2 U_L to the power of-3 [g cm to the power of-1 s to the power of-2]|
| PartType5  |      BH_SoundSpeed                       |Black hole surrounding gas sound speed. Physical c_snd = C_snd U_V [cm/s]|
| PartType5  |      BH_SurroundingGasVel                |Velocity of the gas surrounding the BH (kernel weighted). Physical Velocity = Velocity a to the power of-1 U_M U_V to the power of [cm/s]|
| PartType5  |      BH_TimeLastMerger                   |Expansion factor a when BH particle last accreted an other BH. 0 if the particle as never accreted another BH.                                                       |
| PartType5  |      BH_WeightedDensity                  |Co-moving weighted black hole densities. Physical rho = Densities h to the power of2 a to the power of-3 U_M U_L to the power of-3 [g/cm to the power of3]|
| PartType5  |      Coordinates                         |Co-moving coordinates. Physical: r = a x = Coordinates h to the power of-1 a U_L [cm]|
| PartType5  |      GroupNumber                         |FoF group number particle is in|
| PartType5  |      HostHalo_TVir_Mass                  |Estimate of halo's virial temperature, calculated from the DM halo mass.   T_vir = (MEANMOLIONIZED * PROTONMASS / 3. / BOLTZMANN) * (G * m200 * H(z)) to the power of(2./3.) [K]|
| PartType5  |      Mass                                |Particle mass. Physical m = Mass h to the power of-1 U_M [g]                                                                                                                     |
| PartType5  |      ParticleIDs                         |Unique particle identifier|
| PartType5  |      SmoothingLength                     |Co-moving smoothing length. Physical h = SmoothingLength h to the power of-1 a U_L [cm]|
| PartType5  |      SubGroupNumber                      |Subgroup number particle is in|
| PartType5  |      Velocity                            |Co-moving velocities. Physical v_p = a dx/dt  = Velocities a to the power of1/2 U_V [cm/s]|

### FOF

FOF group numbers begin at 1 and are sorted by mass in descending order. Spherical overdensities associated with a FOF group use the negative of the group number. A particle can only be associated with one spherical overdensity.

The following arrays can be found in the FOF group files:

| Group | Array Name | Description |
| --- | --------------------------------------- | --- |
| FOF | BH_Mdot                                 |Particle mass. Physical m = Mass h to the power of-1 U_M [g] | 
| FOF | BlackHoleMass                           |Particle mass. Physical m = Mass h to the power of-1 U_M [g] | 
| FOF | CentreOfMass                            |Co-moving coordinates. Physical: r = a x = Coordinates h to the power of-1 a U_L [cm] | 
| FOF | GroupLength                             |Number of particles in this group | 
| FOF | GroupLengthType                         |Number of particles in this group of a given type | 
| FOF | GroupMassType                           |Particle mass. Physical m = Mass h to the power of-1 U_M [g] | 
| FOF | GroupOffset                             |Offset of IDs of this group, starts at 0 | 
| FOF | GroupOffsetType                         |Meaning of this variable has not yet been defined. | 
| FOF | Mass                                    |Particle mass. Physical m = Mass h to the power of-1 U_M [g] | 
| FOF | NSF/AExpMaximumTemperature              |Expansion factor a when particle had highest temperature | 
| FOF | NSF/ElementAbundances/Carbon            | | 
| FOF | NSF/ElementAbundances/Helium            | | 
| FOF | NSF/ElementAbundances/Hydrogen          | | 
| FOF | NSF/ElementAbundances/Iron              | | 
| FOF | NSF/ElementAbundances/Magnesium         | | 
| FOF | NSF/ElementAbundances/Neon              | | 
| FOF | NSF/ElementAbundances/Nitrogen          | | 
| FOF | NSF/ElementAbundances/Oxygen            | | 
| FOF | NSF/ElementAbundances/Silicon           | | 
| FOF | NSF/Entropy                             |Meaning of this variable has not yet been defined. | 
| FOF | NSF/Mass                                |Particle mass. Physical m = Mass h to the power of-1 U_M [g] | 
| FOF | NSF/MaximumTemperature                  |Maximum temperature ever reached by a particle [K] | 
| FOF | NSF/Metallicity                         |Mass fraction of elements heavier than Helium | 
| FOF | NSF/SmoothedElementAbundances/Carbon    | | 
| FOF | NSF/SmoothedElementAbundances/Helium    | | 
| FOF | NSF/SmoothedElementAbundances/Hydrogen  | | 
| FOF | NSF/SmoothedElementAbundances/Iron      | | 
| FOF | NSF/SmoothedElementAbundances/Magnesium | | 
| FOF | NSF/SmoothedElementAbundances/Neon      | | 
| FOF | NSF/SmoothedElementAbundances/Nitrogen  | | 
| FOF | NSF/SmoothedElementAbundances/Oxygen    | | 
| FOF | NSF/SmoothedElementAbundances/Silicon   | | 
| FOF | NSF/SmoothedIronMassFracFromSNIa        |Smoothed mass from SNIa divided by particle mass. (Initial particle mass for stars) | 
| FOF | NSF/SmoothedMetallicity                 |Smoothed mass fraction of elements heavier than Helium | 
| FOF | NSF/Temperature                         |Meaning of this variable has not yet been defined. | 
| FOF | ParticleIDs                             |Unique particle identifier | 
| FOF | SF/AExpMaximumTemperature               |Expansion factor a when particle had highest temperature | 
| FOF | SF/ElementAbundaces/Carbon              | | 
| FOF | SF/ElementAbundaces/Helium              | | 
| FOF | SF/ElementAbundaces/Hydrogen            | | 
| FOF | SF/ElementAbundaces/Iron                | | 
| FOF | SF/ElementAbundaces/Magnesium           | | 
| FOF | SF/ElementAbundaces/Neon                | | 
| FOF | SF/ElementAbundaces/Nitrogen            | | 
| FOF | SF/ElementAbundaces/Oxygen              | | 
| FOF | SF/ElementAbundaces/Silicon             | | 
| FOF | SF/Entropy                              |Meaning of this variable has not yet been defined. | 
| FOF | SF/IronMassFracFromSNIa                 |Iron mass from SNIa divided by particle mass. (Initial particle mass for stars) | 
| FOF | SF/Mass                                 |Particle mass. Physical m = Mass h to the power of-1 U_M [g] | 
| FOF | SF/MaximumTemperature                   |Maximum temperature ever reached by a particle [K] | 
| FOF | SF/Metallicity                          |Mass fraction of elements heavier than Helium | 
| FOF | SF/SmoothedElementAbundances/Carbon     | | 
| FOF | SF/SmoothedElementAbundances/Helium     | | 
| FOF | SF/SmoothedElementAbundances/Hydrogen   | | 
| FOF | SF/SmoothedElementAbundances/Iron       | | 
| FOF | SF/SmoothedElementAbundances/Magnesium  | | 
| FOF | SF/SmoothedElementAbundances/Neon       | | 
| FOF | SF/SmoothedElementAbundances/Nitrogen   | | 
| FOF | SF/SmoothedElementAbundances/Oxygen     | | 
| FOF | SF/SmoothedElementAbundances/Silicon    | | 
| FOF | SF/SmoothedIronMassFracFromSNIa         |Smoothed mass from SNIa divided by particle mass. (Initial particle mass for stars) | 
| FOF | SF/SmoothedMetallicity                  |Smoothed mass fraction of elements heavier than Helium | 
| FOF | SF/Temperature                          |Meaning of this variable has not yet been defined. | 
| FOF | StarFormationRate                       |Gas star formation rate in solar masses / yr | 
| FOF | Stars/AExpMaximumTemperature            |Expansion factor a when particle had highest temperature | 
| FOF | Stars/ElementAbundances/Carbon          | | 
| FOF | Stars/ElementAbundances/Helium          | | 
| FOF | Stars/ElementAbundances/Hydrogen        | | 
| FOF | Stars/ElementAbundances/Iron            | | 
| FOF | Stars/ElementAbundances/Magnesium       | | 
| FOF | Stars/ElementAbundances/Neon            | | 
| FOF | Stars/ElementAbundances/Nitrogen        | | 
| FOF | Stars/ElementAbundances/Oxygen          | | 
| FOF | Stars/ElementAbundances/Silicon         | | 
| FOF | Stars/InitialMass                       |Star particle mass at formation time. Physical m = InitialMass h to the power of-1 U_M [g] | 
| FOF | Stars/InitialMassWeightedStellarAge     |Expansion factor a when star particle was born | 
| FOF | Stars/IronMassFracFromSNIa              |Iron mass from SNIa divided by particle mass. (Initial particle mass for stars) | 
| FOF | Stars/Mass                              |Particle mass. Physical m = Mass h to the power of-1 U_M [g] | 
| FOF | Stars/MaximumTemperature                |Maximum temperature ever reached by a particle [K] | 
| FOF | Stars/Metallicity                       |Mass fraction of elements heavier than Helium | 
| FOF | Stars/SmootheElementAbundances/Carbon   | | 
| FOF | Stars/SmootheElementAbundances/Helium   | | 
| FOF | Stars/SmootheElementAbundances/Hydrogen | | 
| FOF | Stars/SmootheElementAbundances/Iron     | | 
| FOF | Stars/SmootheElementAbundances/Magnesium| | 
| FOF | Stars/SmootheElementAbundances/Neon     | | 
| FOF | Stars/SmootheElementAbundances/Nitrogen | | 
| FOF | Stars/SmootheElementAbundances/Oxygen   | | 
| FOF | Stars/SmootheElementAbundances/Silicon  | | 
| FOF | Stars/SmoothedIronMassFracFromSNIa      |Smoothed mass from SNIa divided by particle mass. (Initial particle mass for stars) | 
| FOF | Stars/SmoothedMetallicity               |Smoothed mass fraction of elements heavier than Helium | 
| FOF | Velocity                                |Co-moving velocities. Physical v_p = a dx/dt  = Velocities a to the power of1/2 U_V [cm/s](A bug affects this quantity. The correct factor is a^-1 and not a^0.5) | 

### Subfind

Subgroup numbers begin at 0. Subgroup 0 of a FOF group corresponds to the most massive subgroup within the group.

The following arrays can be found in the Eagle subfind files: 

| Group | Array Name | Description |
| --- | ----------------------------------------------- | --- |
| FOF | ContaminationMass                               |Contaminating mass. Physical M = Mass U_M[g] |  
| FOF | FirstSubhaloID                                  |Index of first sub halo in SubHalo list (starts at 0) |  
| FOF | GroupCentreOfPotential                          |Co-moving position of most bound particle. Physical position = position h to the power of-1 a U_L[cm] |  
| FOF | GroupLength                                     |Number of particles in this group |  
| FOF | GroupMass                                       | Total mass of FoF group. Physical M = Mass h to the power of-1 U_M[g] |  
| FOF | GroupOffset                                     |Offset of IDs of this group, starts at 0 |  
| FOF | Group_M_Crit200                                 | Mass within Rcrit200. Physical M = Mass h to the power of-1 U_M[g] |  
| FOF | Group_M_Crit2500                                |M_Crit2500 |  
| FOF | Group_M_Crit500                                 |M_Crit500 |  
| FOF | Group_M_Mean200                                 |Mass within RMean200. Physical M = Mass h to the power of-1 U_M[g] |  
| FOF | Group_M_Mean2500                                |M_Mean2500 |  
| FOF | Group_M_Mean500                                 |M_Mean500 |  
| FOF | Group_M_TopHat200                               |Mass within RTophat200. Physical M = Mass h to the power of-1 U_M[g] |  
| FOF | Group_R_Crit200                                 |Co-moving radius within which density is 200 times critical density. Physical radius = radius h to the power of-1 a U_L[cm] |  
| FOF | Group_R_Crit2500                                |R_Crit2500 |  
| FOF | Group_R_Crit500                                 |R_Crit500 |  
| FOF | Group_R_Mean200                                 | Co-moving radius within which density is 200 times mean density. Physical radius = radius h to the power of-1 a U_L[cm] |  
| FOF | Group_R_Mean2500                                |R_Mean2500 |  
| FOF | Group_R_Mean500                                 |R_Mean500 |  
| FOF | Group_R_TopHat200                               | Co-moving radius within which density is 200 times (18 * pi to the power of2 + 82 * (Omega_m(z)-1) - 39 * (Omega_m(z)-1) to the power of2); Physical radius = radius h to the power of-1 a U_L[cm] |  
| FOF | NumOfSubhalos                                   |Number of subhaloes in this FoF group |  
| IDs | ParticleID                                      |PID |  
| IDs | Particle_Binding_Energy                         |Binding energy of particles |  
| Subhalo | ApertureMeasurements/Mass/001kpc            |Masses within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/Mass/003kpc            |Masses within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/Mass/005kpc            |Masses within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/Mass/010kpc            |Masses within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/Mass/020kpc            |Masses within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/Mass/030kpc            |Masses within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/Mass/040kpc            |Masses within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/Mass/050kpc            |Masses within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/Mass/070kpc            |Masses within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/Mass/100kpc            |Masses within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/SFR/001kpc             |Star formation rate within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/SFR/003kpc             |Star formation rate within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/SFR/005kpc             |Star formation rate within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/SFR/010kpc             |Star formation rate within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/SFR/020kpc             |Star formation rate within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/SFR/030kpc             |Star formation rate within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/SFR/040kpc             |Star formation rate within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/SFR/050kpc             |Star formation rate within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/SFR/070kpc             |Star formation rate within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/SFR/100kpc             |Star formation rate within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo |  
| Subhalo | ApertureMeasurements/VelDisp/001kpc         |Stellar velocity dispersion within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo ((Expansion factor, aexp-scale-exponent, is incorrect, should be 0.)) |  
| Subhalo | ApertureMeasurements/VelDisp/003kpc         |Stellar velocity dispersion within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo 1) |  
| Subhalo | ApertureMeasurements/VelDisp/005kpc         |Stellar velocity dispersion within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo 1) |  
| Subhalo | ApertureMeasurements/VelDisp/010kpc         |Stellar velocity dispersion within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo 1) |  
| Subhalo | ApertureMeasurements/VelDisp/020kpc         |Stellar velocity dispersion within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo 1) |  
| Subhalo | ApertureMeasurements/VelDisp/030kpc         |Stellar velocity dispersion within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo 1) |  
| Subhalo | ApertureMeasurements/VelDisp/040kpc         |Stellar velocity dispersion within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo 1) |  
| Subhalo | ApertureMeasurements/VelDisp/050kpc         |Stellar velocity dispersion within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo 1) |  
| Subhalo | ApertureMeasurements/VelDisp/070kpc         |Stellar velocity dispersion within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo 1) |  
| Subhalo | ApertureMeasurements/VelDisp/100kpc         |Stellar velocity dispersion within apertures from 1 to 100 kpc (physical) of the centre of potential of the subhalo 1) |  
| Subhalo | BlackHoleMass                               |BH mass. Physical m = Mass h to the power of-1 U_M [g] |  
| Subhalo | BlackHoleMassAccretionRate                  |BH accretion rate. Physical mdot = BH_Mdot h to the power of-1 U_M /U_T [g/s] |  
| Subhalo | CentreOfMass                                |Co-moving position of COM. Physical position = position h to the power of-1 a U_L[cm] |  
| Subhalo | CentreOfPotential                           |Co-moving position of most bound particle. Physical position = position h to the power of-1 a U_L[cm] |  
| Subhalo | GasSpin                                     |Angular momentum per unit mass of gas particles |  
| Subhalo | GroupNumber                                 |FOF Group Number subhalo belongs to |  
| Subhalo | HalfMassProjRad                             |Projected (av. over 3 axes) radius enclosing half of the subhalo mass comprised by each particle type 1) |  
| Subhalo | HalfMassRad                                 |Radius enclosing half of the subhalo mass comprised by each particle type 1) |  
| Subhalo | IDMostBound                                 |Particle ID of lowest *total* energy particle |  
| Subhalo | InertiaTensor                               |Matrix for the second moment of matter distribution. 1) |  
| Subhalo | InitialMassWeightedBirthZ                   |Initial mass weighted metallicity of stars |  
| Subhalo | InitialMassWeightedStellarAge               |Initial mass weighted age of stars in Gyr |  
| Subhalo | KineticEnergy                               |Total kinetic energy of particles bound to this halo |  
| Subhalo | Mass                                        |Total mass of this group. Physical M = Mass h to the power of-1 U_M[g] |  
| Subhalo | MassTwiceHalfMassRad                        |Mass contained within twice the half mass radius for each particle type |  
| Subhalo | MassType                                    |Total mass of this group for each particle type. Physical M = Mass h to the power of-1 U_M[g] |  
| Subhalo | NSF/ElementAbundances/Carbon                |((Some values are incorrect |  
| Subhalo | NSF/ElementAbundances/Helium                | 2) |  
| Subhalo | NSF/ElementAbundances/Hydrogen              | 2) |  
| Subhalo | NSF/ElementAbundances/Iron                  | 2) |  
| Subhalo | NSF/ElementAbundances/Magnesium             | 2) |  
| Subhalo | NSF/ElementAbundances/Neon                  | 2) |  
| Subhalo | NSF/ElementAbundances/Nitrogen              | 2) |  
| Subhalo | NSF/ElementAbundances/Oxygen                | 2) |  
| Subhalo | NSF/ElementAbundances/Silicon               | 2) |  
| Subhalo | NSF/IronFromSNIa                            |Iron from SNIa |  
| Subhalo | NSF/IronFromSNIaSmoothed                    |Smoothed iron from SNIa |  
| Subhalo | NSF/KineticEnergy                           |Kinetic energy of NSF gas |  
| Subhalo | NSF/Mass                                    |Mass |  
| Subhalo | NSF/MassFromAGB                             |Mass from AGB |  
| Subhalo | NSF/MassFromSNII                            |Mass from SNII |  
| Subhalo | NSF/MassFromSNIa                            |Mass from SNIa |  
| Subhalo | NSF/MassWeightedEntropy                     |Mass weighted mean entropy of NSF gas |  
| Subhalo | NSF/MassWeightedPotential                   |Mass weighted potential |
| Subhalo | NSF/MassWeightedTemperature                 |Mass weighted mean temperature of NSF gas |  
| Subhalo | NSF/Metallicity                             |Metallicity weighted by mass of non star forming gas in subhalos |  
| Subhalo | NSF/MetalsFromAGB                           |Mass in metals from AGB |  
| Subhalo | NSF/MetalsFromSNII                          |Mass in metals from SNII |  
| Subhalo | NSF/MetalsFromSNIa                          |Mass in metals from SNIa |  
| Subhalo | NSF/SmootheElementAbuncdances/Carbon        | 2) |  
| Subhalo | NSF/SmootheElementAbuncdances/Helium        | 2) |  
| Subhalo | NSF/SmootheElementAbuncdances/Hydrogen      | 2) |  
| Subhalo | NSF/SmootheElementAbuncdances/Iron          | 2) |  
| Subhalo | NSF/SmootheElementAbuncdances/Magnesium     | 2) |  
| Subhalo | NSF/SmootheElementAbuncdances/Neon          | 2) |  
| Subhalo | NSF/SmootheElementAbuncdances/Nitrogen      | 2) |  
| Subhalo | NSF/SmootheElementAbuncdances/Oxygen        | 2) |  
| Subhalo | NSF/SmootheElementAbuncdances/Silicon       | 2) |  
| Subhalo | NSF/SmoothedMetallicity                     |Smoothed metallicity weighted by mass of non star forming gas in subhalos |  
| Subhalo | NSF/Spin                                    |Angular momentum per unit mass of non-star-forming gas particles |  
| Subhalo | NSF/ThermalEnergy                           |Thermal energy of NSF gas |  
| Subhalo | NSF/TotalEnergy                             |Total energy of NSF gas |  
| Subhalo | SF/ElementAbundances/Carbon                 | 2) |  
| Subhalo | SF/ElementAbundances/Helium                 | 2) |  
| Subhalo | SF/ElementAbundances/Hydrogen               | 2) |  
| Subhalo | SF/ElementAbundances/Iron                   | 2) |  
| Subhalo | SF/ElementAbundances/Magnesium              | 2) |  
| Subhalo | SF/ElementAbundances/Neon                   | 2) |  
| Subhalo | SF/ElementAbundances/Nitrogen               | 2) |  
| Subhalo | SF/ElementAbundances/Oxygen                 | 2) |  
| Subhalo | SF/ElementAbundances/Silicon                | 2) |  
| Subhalo | SF/IronFromSNIa                             |Iron from SNIa |  
| Subhalo | SF/IronFromSNIaSmoothed                     |Smoothed iron from SNIa |  
| Subhalo | SF/KineticEnergy                            |Kinetic energy of SF gas |  
| Subhalo | SF/Mass                                     |Mass |  
| Subhalo | SF/MassFromAGB                              |Mass from AGB |  
| Subhalo | SF/MassFromSNII                             |Mass from SNII |  
| Subhalo | SF/MassFromSNIa                             |Mass from SNIa |  
| Subhalo | SF/MassWeightedEntropy                      |Mass weighted mean entropy of SF gas |  
| Subhalo | SF/MassWeightedPotential                    |Mass weighted potential |  
| Subhalo | SF/MassWeightedTemperature                  |Mass weighted mean temperature of SF gas |  
| Subhalo | SF/Metallicity                              |Metallicity weighted by mass of star forming gas in subhalos |  
| Subhalo | SF/MetalsFromAGB                            |Mass in metals from AGB |  
| Subhalo | SF/MetalsFromSNII                           |Mass in metals from SNII |  
| Subhalo | SF/MetalsFromSNIa                           |Mass in metals from SNIa |  
| Subhalo | SF/SFWeightedMetallicity                    |Metallicity weighted by star formation rate of star forming gas in subhalos |  
| Subhalo | SF/SmoothedElementAbundances/Carbon         | 2) |  
| Subhalo | SF/SmoothedElementAbundances/Helium         | 2) |  
| Subhalo | SF/SmoothedElementAbundances/Hydrogen       | 2) |  
| Subhalo | SF/SmoothedElementAbundances/Iron           | 2) |  
| Subhalo | SF/SmoothedElementAbundances/Magnesium      | 2) |  
| Subhalo | SF/SmoothedElementAbundances/Neon           | 2) |  
| Subhalo | SF/SmoothedElementAbundances/Nitrogen       | 2) |  
| Subhalo | SF/SmoothedElementAbundances/Oxygen         | 2) |  
| Subhalo | SF/SmoothedElementAbundances/Silicon        | 2) |  
| Subhalo | SF/SmoothedMetallicity                      |Smoothed metallicity weighted by mass of star forming gas in subhalos |  
| Subhalo | SF/SmoothedSFWeightedMetallicity            |Smoothed metallicity weighted by star formation rate of star forming gas in subhalos |  
| Subhalo | SF/Spin                                     |Angular momentum per unit mass of star-forming gas particles |  
| Subhalo | SF/ThermalEnergy                            |Thermal energy of SF gas |  
| Subhalo | SF/TotalEnergy                              |Total energy of SF gas |  
| Subhalo | StarFormationRate                           |Total gas star formation rate in solar masses / yr |  
| Subhalo | Stars/ElementAbundances/Carbon              | 2) |  
| Subhalo | Stars/ElementAbundances/Helium              | 2) |  
| Subhalo | Stars/ElementAbundances/Hydrogen            | 2) |  
| Subhalo | Stars/ElementAbundances/Iron                | 2) |  
| Subhalo | Stars/ElementAbundances/Magnesium           | 2) |  
| Subhalo | Stars/ElementAbundances/Neon                | 2) |  
| Subhalo | Stars/ElementAbundances/Nitrogen            | 2) |  
| Subhalo | Stars/ElementAbundances/Oxygen              | 2) |  
| Subhalo | Stars/ElementAbundances/Silicon             | 2) |  
| Subhalo | Stars/IronFromSNIa                          |Iron from SNIa |  
| Subhalo | Stars/IronFromSNIaSmoothed                  |Smoothed iron from SNIa |  
| Subhalo | Stars/KineticEnergy                         |Kinetic energy of stars |  
| Subhalo | Stars/Mass                                  |Mass in stars |  
| Subhalo | Stars/MassFromAGB                           |Mass from AGB |  
| Subhalo | Stars/MassFromSNII                          |Mass from SNII |  
| Subhalo | Stars/MassFromSNIa                          |Mass from SNIa |  
| Subhalo | Stars/MassWeightedPotential                 |Mass weighted potential |  
| Subhalo | Stars/Metallicity                           |Metallicity weighted by mass of stars in subhalos |  
| Subhalo | Stars/MetalsFromAGB                         |Mass in metals from AGB |  
| Subhalo | Stars/MetalsFromSNII                        |Mass in metals from SNII |  
| Subhalo | Stars/MetalsFromSNIa                        |Mass in metals from SNIa |  
| Subhalo | Stars/SmoothedElementAbundances/Carbon      | 2) |  
| Subhalo | Stars/SmoothedElementAbundances/Helium      | 2) |  
| Subhalo | Stars/SmoothedElementAbundances/Hydrogen    | 2) |  
| Subhalo | Stars/SmoothedElementAbundances/Iron        | 2) |  
| Subhalo | Stars/SmoothedElementAbundances/Magnesium   | 2) |  
| Subhalo | Stars/SmoothedElementAbundances/Neon        | 2) |  
| Subhalo | Stars/SmoothedElementAbundances/Nitrogen    | 2) |  
| Subhalo | Stars/SmoothedElementAbundances/Oxygen      | 2) |  
| Subhalo | Stars/SmoothedElementAbundances/Silicon     | 2) |  
| Subhalo | Stars/SmoothedMetallicity                   |Smoothed metallicity weighted by mass of stars in subhalos |  
| Subhalo | Stars/Spin                                  |Angular momentum per unit mass of star particles |  
| Subhalo | Stars/TotalEnergy                           |Total energy of stars |  
| Subhalo | StellarInitialMass                          |Stellar initial mass |  
| Subhalo | StellarVelDisp                              |Stellar velocity dispersion <sup>1)</sup>|  
| Subhalo | StellarVelDisp_HalfMassProjRad              |Stellar velocity dispersion within half mass radius <sup>1)</sup>|  
| Subhalo | SubGroupNumber                              |SubGroup Number of subhalo, begins at 0 for most massive subhalo within a group |  
| Subhalo | SubLength                                   |Number of particles in this subhalo |  
| Subhalo | SubLengthType                               |Number of particles of each type in this subhalo |  
| Subhalo | SubOffset                                   |Offset of IDs in this subhalo.  Starts at 0 |  
| Subhalo | ThermalEnergy                               |Total thermal energy of particles bound to this halo |  
| Subhalo | TotalEnergy                                 |Total energy of particles bound to this halo |  
| Subhalo | Velocity                                    |Vel |  
| Subhalo | Vmax                                        |Co-moving maximum circular velocity. Physical velocity = ??? |  
| Subhalo | VmaxRadius                                  |VmaxRad |  
| | | |
| footnote | | 1) Expansion factor, aexp-scale-exponent, is incorrect, should be 0. |
| footnote | | 2) Some values are incorrect |
