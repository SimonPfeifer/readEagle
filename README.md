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

 