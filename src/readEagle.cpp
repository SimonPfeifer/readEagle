#include <Python.h>
#include <numpy/ndarrayobject.h>
//#include <numpy/ndarraytypes.h>
#include "numpy/arrayobject.h"
#include "numpy/npy_common.h"
 
#include <string>
#include <iostream>
#include <cctype>

#undef H5T_NATIVE_FLOAT_g

#include "readArray.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"


inline fileType getFileType(const std::string& type)
{
  if(type == "FOF")
    return FOF;
  else if(type == "FOF_PARTICLES")
    return FOF_PARTICLES;
  else if(type == "SNIP_FOF")
    return SNIP_FOF;
  else if(type == "SNIP_FOF_PARTICLES")
    return SNIP_FOF_PARTICLES;
  else if(type == "SNAP")
    return SNAP;
  else if(type == "SNAPSHOT")
    return SNAP;
  else if(type == "SNIP")
    return SNIP;
  else if(type == "SNIPSHOT")
    return SNIP;
  else if(type == "SUBFIND")
    return SUBFIND;
  else if(type == "SUBFIND_GROUP")
    return SUBFIND_GROUP;
  else if(type == "SUBFIND_PARTICLES")
    return SUBFIND_IDS;
  else if(type == "SUBFIND_IDS")
    return SUBFIND_IDS;
  else if(type == "SNIP_SUBFIND")
    return SUBFIND;
  else if(type == "SNIP_SUBFIND_GROUP")
    return SUBFIND_GROUP;
  else if(type == "SNIP_SUBFIND_PARTICLES")
    return SUBFIND_IDS;
  else if(type == "SNIP_SUBFIND_IDS")
    return SUBFIND_IDS;
  else if(type == "PARTDATA")
    return PARTDATA;
  else if(type == "SNIP_PARTDATA")
    return PARTDATA;
  else
    error("Invalid fileType provided.");
} 

template<typename T>inline char* pyFormatter();
template<>inline char* pyFormatter<int>(){ return "i";}
template<>inline char* pyFormatter<unsigned int>(){ return "I";}
template<>inline char* pyFormatter<unsigned long long>(){ return "K";}
template<>inline char* pyFormatter<float>(){ return "f";}
template<>inline char* pyFormatter<double>(){ return "d";}

template<typename T> inline NPY_TYPES pyType();
template<>inline NPY_TYPES pyType<int>(){ return NPY_INT;}
template<>inline NPY_TYPES pyType<unsigned int>(){ return NPY_UINT;}
template<>inline NPY_TYPES pyType<unsigned long long>(){ return NPY_ULONGLONG;}
template<>inline NPY_TYPES pyType<float>(){ return NPY_FLOAT;}
template<>inline NPY_TYPES pyType<double>(){ return NPY_DOUBLE;}
 
template<int size, typename T>
inline PyObject* returnArrayAttribute(fileType ftype, const std::string& dir, const std::string& tag, const std::string& name)
{
  npy_intp dims[1] = {size};
  import_array()
  PyArrayObject* array = (PyArrayObject *) PyArray_SimpleNew(1, dims,  pyType<T>());
  const Tuple<size, T> ret = readAttribute<Tuple<size, T> >(ftype, dir, tag, name);
  for(int i(0); i<size; ++i)
    ((T*) (array->data))[i] = ret[i];
  return (PyObject*) array;
}

template<typename T>
PyObject* returnAttribute(fileType ftype, const std::string& dir, const std::string& tag, const std::string& name, int size)
{
  if(size == 1)
    {
      const T ret = readAttribute<T>(ftype, dir, tag, name);
      return Py_BuildValue(pyFormatter<T>(), ret);
    }
  else
    {
      switch(size)
      	{
      	case 6: return returnArrayAttribute<6,T>(ftype, dir, tag, name);
	default:
	  error("Atrribute size not supported. Buying Matthieu a drink might help.");
      	}
    }
  return 0;
}


PyObject* readAttributePython(PyObject* self, PyObject *args, PyObject *kwargs)
{
  char *c_type, *c_dir, *c_tag, *c_name;
  static char* list[] = {"fileType", "directory", "tag", "attribute", 0};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "ssss", list, &c_type, &c_dir, &c_tag, &c_name ))
    return 0;

  /* Convert to C++ string types for convenience */
  const std::string type=c_type;
  const std::string dir=c_dir;
  const std::string tag=c_tag;
  const std::string name=c_name;

  /* Convert the string to the internal enum value */
  const fileType ftype = getFileType(type);

  /* Recover the type of the attribute */
  const hid_t htype = readAttributeType(ftype, dir, tag, name);

  /* Recover the type of the attribute */
  const int size = readAttributeSize(ftype, dir, tag, name);

  /* Call the version of the function that corresponds to the type */
  if(H5Tequal(htype, H5T_NATIVE_INT))
    {
      H5Tclose(htype);
      return returnAttribute<int>(ftype, dir, tag, name, size);
    }
  else if(H5Tequal(htype, H5T_NATIVE_UINT32))
    {
      H5Tclose(htype);
      return returnAttribute<unsigned int>(ftype, dir, tag, name, size);
    }
  else if(H5Tequal(htype, H5T_NATIVE_UINT64))
    {
      H5Tclose(htype);
      return returnAttribute<unsigned long long>(ftype, dir, tag, name, size);
    }
  else if(H5Tequal(htype, H5T_NATIVE_FLOAT))
    {
      H5Tclose(htype);
      return returnAttribute<float>(ftype, dir, tag, name, size);
    }
  else if(H5Tequal(htype, H5T_NATIVE_DOUBLE))
    {
      H5Tclose(htype);
      return returnAttribute<double>(ftype, dir, tag, name, size);
    }
  else
    {
      error("Unsupported data type in HDF5 file. Buying Matthieu a drink might help.");
    }

  return 0;
}



// --------------------------------------------------------------------------------------------------------------


template<typename T>
PyObject* returnArray(fileType ftype, const std::string& dir, const std::string& tag, const std::string& name, bool physicalUnits, bool noH, bool useCGS, bool verbose, int numThreads, bool oldSubfindFix=false, bool owlsStyle=false)
{
 /* Open file #0 to read global information */  
  const hid_t h_file = H5Fopen(fileName(ftype, dir, tag).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if(h_file < 0)
    error("Impossible to read file #0");

  /* Decide which of the 6 values for the number of elements to read if we are dealing with a sn[ai]pshot */
  const int partType = getGadgetType(ftype, name);

  /* Decide which of the 6 values for the number of elements to read if we are dealing with a sn[ai]pshot */
  const int numSnap = readAttribute<int>(h_file, numFilesFieldName(ftype));
  const size_t totNumElem = readAttribute<Tuple<6, size_t> >(h_file, oldSubfindFix ? "/Header/NumPart_Total" : totElementFieldName(ftype))[oldSubfindFix ? 2 : partType];
  H5Fclose(h_file);

  const int length = readArraySize(ftype, dir, tag, name);

  if(verbose)
    std::cout << "Extracting " << totNumElem << " values for '"<< name << "' (PartType=" << partType << ") in " << numSnap << " " << fileTypeName(ftype) << " files with tag '" << tag << "' using " << numThreads << " thread(s)" << ( oldSubfindFix ? " (with old subfind fix)." : "." ) << std::endl;

  /* Create array */
  import_array()
  PyArrayObject* array = 0;

  if(length == 1)
    {
      npy_intp dims[1] = {static_cast<npy_intp>(totNumElem)};
      array = (PyArrayObject *) PyArray_SimpleNew(1, dims,  pyType<T>());
    }
  else
    {
      npy_intp dims[2] = {static_cast<npy_intp>(totNumElem), length};
      array = (PyArrayObject *) PyArray_SimpleNew(2, dims,  pyType<T>());
    }

  import_array()
  T* c_array = (T*) array->data;

  /* Separate arrays for parallel reads */
  std::vector<size_t> sizes(numSnap, 0);
  std::vector<size_t> offsets(numSnap, 0);

  /* Prepare data for threads*/
  pthread_t* threads = new pthread_t[numThreads];
  thread_info<T>* thread_out = new thread_info<T>[numThreads];
  for(int i(0); i<numThreads; ++i)
    {
      thread_out[i].id = i;
      thread_out[i].numSnap = numSnap;
      thread_out[i].numThreads = numThreads;
      thread_out[i].type = ftype;
      thread_out[i].partType = partType;
      thread_out[i].dir = dir;
      thread_out[i].tag = tag;
      thread_out[i].arrayName = name;
      thread_out[i].array = &c_array[0];
      thread_out[i].sizes = &sizes[0];
      thread_out[i].offsets = &offsets[0];
      thread_out[i].oldSubfindFix = oldSubfindFix;
    }

  timeval begin, end;
  gettimeofday(&begin, NULL);

  if(verbose)
    std::cout << "Reading array lengths" << std::flush;
  
 
  /* Compute sizes for parallel reads using threads */
  for(int i(0); i<numThreads; ++i)
    {
      const int ret = pthread_create(&threads[i], 0, thread_readSize<T>, &thread_out[i]);
      if(ret != 0)
	error("Error while launching thread");
    }

  for (int i(0); i < numThreads; ++i) 
    {
      const int ret = pthread_join(threads[i], 0);
      if(ret != 0)
	error("Error while stopping thread");
    } 
 
  gettimeofday(&end, NULL);  
  if(verbose)
    std::cout << " took " << ((end.tv_sec-begin.tv_sec)*1000000 + end.tv_usec-begin.tv_usec) / 1e6 << "s" << std::endl;

  // for(int i(0); i<numSnap; ++i)
  //   sizes[i] *= length;

  /* Do a scan to convert from size to offset */
  for(int i(1); i<numSnap; ++i)
      offsets[i] = offsets[i-1] + sizes[i-1];

  if( offsets[numSnap-1] + sizes[numSnap-1] != totNumElem )
    error( "Sizes mismatch !" );
 
  /* Read data from files using threads */
  gettimeofday(&begin, NULL);
  if(verbose) 
    std::cout << "Reading data" << std::flush;

  for(int i(0); i<numThreads; ++i)
    {
      int ret = -1;
      switch(length)
	{
	case 1:
	  ret = pthread_create(&threads[i], 0, thread_readData<T, T >, &thread_out[i]); break;
	case 3:
	  ret = pthread_create(&threads[i], 0, thread_readData<T, Vector<T> >, &thread_out[i]); break;
	case 6:
	  ret = pthread_create(&threads[i], 0, thread_readData<T, Tuple<6,T> >, &thread_out[i]); break;
	default:
	  error("Array size not supported. Talk to Matthieu...");
	}
      if(ret != 0)
	error("Error while launching thread");
    }

  for (int i(0); i < numThreads; ++i) 
    {
      const int ret = pthread_join(threads[i], 0);
      if(ret != 0)
	error("Error while stopping thread");
    } 

  gettimeofday(&end, NULL);  
  if(verbose)
    std::cout << " took " << ((end.tv_sec-begin.tv_sec)*1000000 + end.tv_usec-begin.tv_usec) / 1e6 << "s" << std::endl;

  /* Find a non-empty file for unit conversion */
  int nonEmpty = -1;
  for(int i(0); i<numSnap; ++i)
    {
      if(sizes[i] > 0)
	{
	  nonEmpty = i;
	  break;
	}
    }

  if(physicalUnits || noH || useCGS)
    {

      /* Read factors for unit conversion */
      const hid_t h_file = H5Fopen(fileName(ftype, dir, tag).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      const double hubbleParam = readAttribute<double>(h_file, "/Header/HubbleParam");
      const double a = readAttribute<double>(h_file, "/Header/ExpansionFactor");
      H5Fclose(h_file);
      const hid_t h_file2 = H5Fopen(fileName(ftype, dir, tag, nonEmpty).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      const float a_exponant = readUnits<float>(h_file2, name + "/aexp-scale-exponent");
      const float h_exponant = readUnits<float>(h_file2, name + "/h-scale-exponent");
      const double cgsConversion = readUnits<double>(h_file2, name + "/CGSConversionFactor");
      H5Fclose(h_file2);
      
      /* Apply comoving to physical units conversion */
      if(physicalUnits)
	{
	  if(a_exponant != 0.)
	    {
	      if(verbose)
		std::cout << "Converting to physical units. (Multiplication by a^" << a_exponant << ", a=" << a  << ")" << std::endl;
	      const double factor = std::pow(a, static_cast<double>(a_exponant));
	      for(size_t i(0); i<totNumElem*length; ++i)
		c_array[i] *= factor;
	    }
	  else
	    {
	      if(verbose)
		std::cout << "Converting to physical units. No conversion needed !" << std::endl;
	    }
	}
      
      /* Apply h units conversion */
      if(noH)
	{
	  if(h_exponant != 0.)
	    {
	      if(verbose)
		std::cout << "Converting to h-free units. (Multiplication by h^" << h_exponant << ", h=" << hubbleParam << ")" << std::endl;
	      const double factor = std::pow(hubbleParam, static_cast<double>(h_exponant));
	      for(size_t i(0); i<totNumElem*length; ++i)
		c_array[i] *= factor;
	    }
	  else
	    {
	      if(verbose)
		std::cout << "Converting to h-free units. No conversion needed !" << std::endl;
	    }
	}
      
      /* Apply CGS units conversion */
      if(useCGS)
	{
	  if(cgsConversion != 0.)
	    {
	      if(verbose)
		std::cout << "Converting to CGS units. (Multiplication by " << cgsConversion << ")" << std::endl;
	      for(size_t i(0); i<totNumElem*length; ++i)
		c_array[i] *= cgsConversion;
	    }
	  else
	    {
	      if(verbose)
		std::cout << "Converting to CGS units. No conversion needed !" << std::endl;
	    }
	}
    }


  delete[] threads;
  delete[] thread_out;


  return (PyObject*)array;
}


PyObject* readArrayPython(PyObject* self, PyObject *args, PyObject *kwargs)
{
  char *c_type, *c_dir, *c_tag, *c_name;
  bool physicalUnits = true , noH = true, useCGS = false, verbose = true, oldSubfindFix = false, owlsStyle = false;
  int numThreads = 16;
  static char* list[] = {"fileType", "directory", "tag", "attribute", "physicalUnits", "noH", "useCGS", "verbose", "numThreads", "oldSubfindFix", "OWLSstyle", 0};
  if(!PyArg_ParseTupleAndKeywords(args, kwargs, "ssss|iiiiiii", list, &c_type, &c_dir, &c_tag, &c_name, &physicalUnits, &noH, &useCGS, &verbose, &numThreads, &oldSubfindFix, &owlsStyle ))
    return 0;

  /* Convert to C++ string types for convenience */
  const std::string type=c_type;
  const std::string dir=c_dir;
  const std::string tag=c_tag;
  const std::string name=c_name;

  /* Check that we have at least one thread */
  if(numThreads < 1)
    error("Invalid number of threads. Value must be >0.");

  /* Convert the string to the internal enum value */
  const fileType ftype = getFileType(type);

  /* Recover the type of the arrray */
  const hid_t htype = readArrayType(ftype, dir, tag, name);

  /* Check validity of fix */
  if(oldSubfindFix && ftype != SUBFIND_IDS)
    error("Can only use subfind fix when reading subfind IDs.");

  /* Call the version of the function that corresponds to the type */
  if(H5Tequal(htype, H5T_NATIVE_INT))
    {
      H5Tclose(htype);
      return returnArray<int>(ftype, dir, tag, name, physicalUnits, noH, useCGS, verbose, numThreads, oldSubfindFix, owlsStyle);
    }
  else if(H5Tequal(htype, H5T_NATIVE_UINT32))
    {
      H5Tclose(htype);
      return returnArray<unsigned int>(ftype, dir, tag, name, physicalUnits, noH, useCGS, verbose, numThreads, oldSubfindFix, owlsStyle);
    }
  else if(H5Tequal(htype, H5T_NATIVE_UINT64))
    {
      H5Tclose(htype);
      return returnArray<unsigned long long>(ftype, dir, tag, name, physicalUnits, noH, useCGS, verbose, numThreads, oldSubfindFix, owlsStyle);
    }
  else if(H5Tequal(htype, H5T_NATIVE_FLOAT))
    {
      H5Tclose(htype);
      return returnArray<float>(ftype, dir, tag, name, physicalUnits, noH, useCGS, verbose, numThreads, oldSubfindFix, owlsStyle);
    }
  else if(H5Tequal(htype, H5T_NATIVE_DOUBLE))
    {
      H5Tclose(htype);
      return returnArray<double>(ftype, dir, tag, name, physicalUnits, noH, useCGS, verbose, numThreads, oldSubfindFix, owlsStyle);
    }
  else
    {
      error("Unsupported data type in HDF5 file. Buying Matthieu a drink might help.");
    }


  return 0;
}

// --------------------------------------------------------------------------------------------------------------


// PyObject* example(PyObject* self, PyObject *args, PyObject *kwargs)
// {
//   int argument1, argument2;
//   float argument3;
//   static char* list[] = {"argument1", "argument2", "argument3", 0};
//   if(!PyArg_ParseTupleAndKeywords(args, kwargs, "iif", list, &argument1, &argument2, &argument3))
//     return 0;

//   std::cout << "my first function !" << std::endl;

//   npy_intp dims[2] = {3,3};
//   PyArrayObject* array = (PyArrayObject *) PyArray_SimpleNew(2, dims,  NPY_DOUBLE);
//   double* c_array = (double*) array->data;

//   int i;
//   for(i=0; i<9; ++i)
//     c_array[i] = i+1;

//   return (PyObject*) array;
// }


static PyMethodDef functions[]={
  //  {"example", (PyCFunction)example, METH_VARARGS|METH_KEYWORDS, "One comment"},
  {"readAttribute", (PyCFunction)readAttributePython, METH_VARARGS|METH_KEYWORDS, "Read one attribute from an EAGLE output file."},
  {"readArray", (PyCFunction)readArrayPython, METH_VARARGS|METH_KEYWORDS, "Read one array from an EAGLE output file."},
  {0, 0, 0, 0}
};

static struct PyModuleDef eagle =
{
    PyModuleDef_HEAD_INIT,
    "eagle", /* name of module */
    "See website - http://eagle.strw.leidenuniv.nl/wiki/doku.php?id=eagle:documentation:reading_python", /* module documentation, may be NULL */
    0, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    functions
};

PyMODINIT_FUNC PyInit_eagle(void)
{
    return PyModule_Create(&eagle);
}

