# readEagle

This python module can be used to read in data from the EALGE and BAHAMAS simultion HDF5 output files. It also applies h-factor corrections and can return the quantities in physical or comoving coordinates. The reading of the data uses an underlying set of C++ routines that read the different files that make up a snapshot in parrallel for higher throughput.

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
 
attr = E.readAttribute(file_type,  directory,  tag,   attribute)
arr = E.readArray(file_type,  directory,  tag,   array)
```

Both routines are constructed in the same way and need to be supplied with the same 4 arguments. The first one, `file_type`, is a string describing the type of file and data read. An example of some allowed arguments are:
```python
file_type = "SNAPSHOT"
file_type = "FOF"
file_type = "SUBFIND"
```

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
attribute = "/Header/Redshift"
array = "/PartType0/Coordinates"
```

The routine returns a numpy array containing the values extracted from the files. The order of the elements is preserved and the type of the values is the same than stored in the HDF5 files. 

There are lists of acceptabe `file_type` and `array` arguments in the `examples/` directory.

### Unit conversion

By default, the `readArray` function converts the data read from the file into “h-free” physical units. This is done by reading the relevant conversion factors from the HDF5 file. The conversions applied to the data are reported by the function and printed to the standard output. For instance, reading the particle coordinates at redshift 1 with this code:
```python
pos = E.readArray("PARTDATA", directory, "019_z001p004", "/PartType0/Coordinates")
```
will yield the following output to stdout:
```
Extracting 6117616 values for '/PartType0/Coordinates' (PartType=0) in 16 ParticleData files with tag '019_z001p004' using 16 thread(s).
Reading array lengths took 0.35501s
Reading data took 3.83977s
Converting to physical units. (Multiplication by a^1, a=0.498972)
Converting to h-free units. (Multiplication by h^-1, h=0.6777)
```

This implies that the data in the file (comoving h-full units) has been converted to physical h-free units by the function. No further conversion is required. This relies on the fact that the units written in the file are correct. **Always check that this is the case by looking at the standard output!**

This behaviour can be modified using the two optional parameters `noH` and `physicalUnits`. If `noH` is set to False then the routine does not apply any h-factor correction. If `physicalUnits` is set to False then no a-factor correction is applied. Running with `noH=False, physicalUnits=False` will hence read in the data as it is in the file without applying any correction. For instance, reading the particle coordinates at redshift 1 with this code:
```python
pos = E.readArray("PARTDATA", directory, "019_z001p004", "/PartType0/Coordinates", noH=False, physicalUnits=False)
```
will yield:
```
Extracting 6117616 values for '/PartType0/Coordinates' (PartType=0) in 16 ParticleData files with tag '019_z001p004' using 16 thread(s).
Reading array lengths took 1.26781s
Reading data took 14.4998s
```
Not displaying (and applying) any correction. The h and a factor corrections can obviously also be applied independently.

### CGS unit conversion

The `readArray` function can also directly convert the data to CGS units. This is done by adding `useCGS=True` to the parameter of the function when calling it. This is independant of the a-factor and h-factor corrections.

Reading the “R_500” radii of all halos in “physical CGS h-free units” would be done using the call:
```python
R_500 = eagle.readArray("SUBFIND_GROUP", directory, tag, "FOF/Group_R_Crit500", useCGS=True)
```
which would lead to the output:
```
Extracting 1885062 values for 'FOF/Group_R_Crit500' (PartType=0) in 256 Subfind (group tabs) files with tag '028_z000p000' using 64 thread(s).
Reading array lengths took 13.646s
Reading data took 42.9861s
Converting to physical units. (Multiplication by a^1, a=1)
Converting to h-free units. (Multiplication by h^-1, h=0.6777)
Converting to CGS units. (Multiplication by 3.08568e+24)
```

### Additional options

There are three additional options to control the behaviour of the `readArray` routine.

The first one, `verbose`, controls whether the function writes what it does to the standard output. Setting `verbose=False` will switch off all the output. Example:
```python
pos = E.readArray("PARTDATA", directory, "019_z001p004", "/PartType0/Coordinates", verbose=False)
```

The second one controls the number of threads used to read the file in parallel. This is done by giving a value to the parameter `numThreads` where the default is 16. Note that increasing the number of threads does not seem to lead to a significant improvement in speed.
```python
pos = E.readArray("PARTDATA", directory, "019_z001p004", "/PartType0/Coordinates") # Read the data using 16 threads
pos = E.readArray("PARTDATA", directory, "019_z001p004", "/PartType0/Coordinates", numThreads=32) # Read the data using 32 threads
```
This only modifies the internal behaviour of the function. Once the function returns, you are given a plain old numpy array without any fancy features and can use it as any other numpy array in a single-threaded python context.

The last one, `oldSubfindFix`, is a workaround for a bug in the older versions of SUBFIND where the number of particle IDs stored was written as an int instead of a long long. Switching this option on allows to read the particle IDs out of old SUBFIND data sets with more than 10^9 particles. All the other quantities can safely be read without this fix.
```python
part_ids = E.readArray("SUBFIND_PARTICLES", directory, "019_z001p004", "IDs/ParticleID", oldSubfindFix=True)
```

There are some examples in the `examples/` directory that show the application of the package as well as lists of acceptable `file_type` and `array` arguments.
