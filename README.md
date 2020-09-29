# ReadEagle

This python module is aimed at reading in data from all the possible EALGE and BAHAMAS hdf5 output files. It also applies h-factor corrections and can return the quantities in physical or comoving coordinates. The reading of the data uses an underlying set of C++ routines that read the different files that make up a snapshot in parrallel for higher throughput.

Disclaimer: There is no reason to believe that the routine presented here will work, return correct results, compile or not create a black hole absorbing the Universe. Please check what you do and verify that the results returned are sensible.

The source code and examples were written by Matthieu Schaller. The python wrapper was converted to Python 3 by Simon Pfeifer. If you have a question, a feature request or a bug report to submit, please create an Issue/Pull request on GitHub or contact spfeifer@aip.de directly.

## Installation

This installation guide is specific to the LJMU HPC computing clusters.

To install the `ReadEagle` package, `ssh` into one of the computing clusters (e.g. famine) and load the HDF5 library using:

```
module load hdf5/1.8.XX
```

I have only tested the installation with HDF5 1.8, where the last version I used was 1.8.21. Using the HDF5 library on the computing clusters rather than a privately installed version has one caveat (saftey feature), and that it that you can only use the `ReadEagle` if you are on one of the clusters. It won't work on "external" or your "star-pc" since it won't be able to find the HDF5 library. This has almost no draw-backs since loading simulation data takes a lot of memory and safes you accidentally crashing you pc.

Once the HDF5 library is loaded, run the following commands:
```
cd src/
python setup.py install --user
```

This will install the module into your personal python directory and assumes that `python` calls a version of Python 3. If more than one version is installed, you can check which version is called by typing `python --version`.

To quickly check if the installation was successful, run `python` in your terminal and import the package using `import eagle`.

## Reading data

The module offers two functions. One to read attributes and one to read arrays:

```python
import eagle as E
 
z = E.readAttribute(fileType,  directory,  tag,   attribute)
M_200 = E.readArray(fileType,  directory,  tag,   array)
```

Both routines are constructed in the same way and need to be supplied with the same 4 arguments. The first one, `fileType`, is a string describing the type of file and data read. The allowed values of this parameter are:  

|fileType | Description | Example of data that can be read |
| --- | --- | --- |
|“FOF” | FoF group informations | Group centre of mass, group length, group star formation rate |
|“FOF_PARTICLES” | IDs of the particles in a FOF group | Particle IDs |
|“SNIP_FOF” | FoF group informations (snipshot) | Group centre of mass, group length, group star formation rate |
|“SNIP_FOF_PARTICLES” | IDs of the particles in a FOF group (snipshot) | Particle IDs |
|“PARTDATA” | Particles that are in a FOF group | Particle Mass, velocity, entropy, stellar age |
|“SNIP_PARTDATA” | Particles that are in a FOF group (snipshot) | Particle Mass, velocity, entropy, stellar age |
|“SNAPSHOT” | Full information about all particles | Particle Mass, velocity, entropy, stellar age |
|“SNIPSHOT” | Reduced information about all particles | Mass, velocity |
|“SUBFIND” | Subhalo information | Subhalo mass, subhalo centre of potential |
|“SUBFIND_GROUP” | Subfind halo information | Group centre of potential, M_200, R_500 |
|“SUBFIND_PARTICLES” | IDs of the particles in a subhalo | Particle IDs |
|“SNIP_SUBFIND” | Subhalo information (snipshot) | Subhalo mass, subhalo centre of potential |
|“SNIP_SUBFIND_GROUP” | Subfind halo information (snipshot) | Group centre of potential, M_200, R_500 |
|“SNIP_SUBFIND_PARTICLES” | IDs of the particles in a subhalo (snipshot) | Particle IDs |


The second argument is the location of the directory containing the data. For instance:
```python
directory = "hpcdata0/simulations/BAHAMAS/RUNNING/HYDRO/L400N1024/nrun/run_m0.00025/data/"
```
Note that not all simulation directories are structured the same way. For example, the `MASSIVE_NEUTRINOS` directory splits the particle data and SUBFIND output into separate folders. That means that you would have to modify the directory depending on whether you are reading data from a snapshot (e.g. `PARTDATA`) or SUBFIND output (e.g. `SUBFIND`), i.e. which `fileType` argument you are using.

The third argument is the `tag` of the output. This is the part of the filename that contains the snapshot number and sometimes the redshift:
```python
tag_without_redshift = "032"
tag_with_redshift = "032_z000p000"
```
Check the simulation directory to see which format is used.


The last argument is the name of the `array` or `attribute` in the HDF5 file structure you want to read. For instance:

```python
attribute = "/Header/BoxSize"
array = "/PartType1/Coordinates"
```

The routine returns a numpy array containing the values extracted from the files. The order of the elements is preserved and the type of the values is the same than stored in the HDF5 files. 