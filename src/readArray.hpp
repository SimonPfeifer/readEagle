#ifndef READ_ARRAY_HPP
#define READ_ARRAY_HPP


#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <cmath>
#include <cassert>
#include <hdf5.h>
#include <sys/time.h>
#include <pthread.h>

#include "vector.hpp"
#include "tuple.hpp"
#include "error.hpp"

enum fileType{FOF,            //FOF group arrays
	      FOF_PARTICLES,  //Particle info in FOF files
	      SNIP_FOF,            //FOF group arrays
	      SNIP_FOF_PARTICLES,  //Particle info in FOF files
	      SNAP, 
	      SNIP, 
	      SUBFIND,        // Subfind sub-halo arrays
	      SUBFIND_GROUP,  // Subfind halo arrays
	      SUBFIND_IDS,
	      SNIP_SUBFIND,        // Subfind sub-halo arrays
	      SNIP_SUBFIND_GROUP,  // Subfind halo arrays
	      SNIP_SUBFIND_IDS,
	      PARTDATA,
	      SNIP_PARTDATA
};

template<typename T> inline hid_t hdf5Type();
template<> inline hid_t hdf5Type<char*>(){return H5T_C_S1; }
template<> inline hid_t hdf5Type<int>(){return H5T_NATIVE_INT; }
template<> inline hid_t hdf5Type<size_t>(){return H5T_NATIVE_UINT64; }
template<> inline hid_t hdf5Type<unsigned int>(){return H5T_NATIVE_UINT32; }
template<> inline hid_t hdf5Type<unsigned long long>(){return H5T_NATIVE_UINT64; }
template<> inline hid_t hdf5Type<float>(){return H5T_NATIVE_FLOAT; }
template<> inline hid_t hdf5Type<double>(){return H5T_NATIVE_DOUBLE; }
template<> inline hid_t hdf5Type<Tuple<6, int> >(){return H5T_NATIVE_INT; }
template<> inline hid_t hdf5Type<Tuple<6, unsigned int> >(){return H5T_NATIVE_UINT32; }
template<> inline hid_t hdf5Type<Tuple<6, unsigned long long> >(){return H5T_NATIVE_UINT64; }
template<> inline hid_t hdf5Type<Tuple<6, size_t> >(){return H5T_NATIVE_UINT64; }
template<> inline hid_t hdf5Type<Tuple<6, float> >(){return H5T_NATIVE_FLOAT; }
template<> inline hid_t hdf5Type<Tuple<6, double> >(){return H5T_NATIVE_DOUBLE; }

/**
 * \brief Returns the name of the file for a given directory, tag, fileNumber and type of file
 */
inline std::string fileName(fileType type, const std::string& dir, const std::string& tag, int i=0)
{
  std::ostringstream oss;
  switch(type)
    {
    case FOF:
    case FOF_PARTICLES:
      oss << dir << "/groups_" << tag << "/group_tab_" << tag << "." << i << ".hdf5";
      return oss.str();
    case SNIP_FOF:
    case SNIP_FOF_PARTICLES:
      oss << dir << "/groups_snip_" << tag << "/group_snip_tab_" << tag << "." << i << ".hdf5";
      return oss.str();
    case SUBFIND:
    case SUBFIND_GROUP:
    case SUBFIND_IDS:
      oss << dir << "/groups_" << tag << "/eagle_subfind_tab_" << tag << "." << i << ".hdf5";
      return oss.str();
    case SNIP_SUBFIND:
    case SNIP_SUBFIND_GROUP:
    case SNIP_SUBFIND_IDS:
      oss << dir << "/groups_snip_" << tag << "/eagle_subfind_snip_tab_" << tag << "." << i << ".hdf5";
      return oss.str();
    case SNIP:
      oss << dir << "/snipshot_" << tag << "/snip_" << tag << "." << i << ".hdf5";
      return oss.str();
    case SNAP:
      oss << dir << "/snapshot_" << tag << "/snap_" << tag << "." << i << ".hdf5";
      return oss.str();
    case PARTDATA:
      oss << dir << "/particledata_" << tag << "/eagle_subfind_particles_" << tag << "." << i << ".hdf5";
      return oss.str();
    case SNIP_PARTDATA:
      oss << dir << "/particledata_snip_" << tag << "/eagle_subfind_snip_particles_" << tag << "." << i << ".hdf5";
      return oss.str();
    default:
      error("Type of files not supported");
    }
}

inline std::string numFilesFieldName(fileType type)
{
  switch(type)
    {
    case FOF:
    case FOF_PARTICLES:
    case SNIP_FOF:
    case SNIP_FOF_PARTICLES:
    case SUBFIND:
    case SUBFIND_GROUP:
    case SUBFIND_IDS:
    case SNIP_SUBFIND:
    case SNIP_SUBFIND_GROUP:
    case SNIP_SUBFIND_IDS:
      return "/Header/NTask";
    case SNAP:
    case SNIP:
    case PARTDATA:
    case SNIP_PARTDATA:
      return "/Header/NumFilesPerSnapshot";
    default:
      error("Type of file unknown");
    }
}


/**
 * \brief Returns the header field containg the total number of elements
 */
inline std::string totElementFieldName(fileType type)
{
  switch(type)
    {
    case FOF_PARTICLES: 
    case SNIP_FOF_PARTICLES: 
      return "/Header/ToTNids";
    case FOF: 
    case SUBFIND_GROUP:
    case SNIP_FOF: 
    case SNIP_SUBFIND_GROUP:
      return "/Header/TotNgroups";
    case SUBFIND: 
    case SNIP_SUBFIND: 
      return "/Header/TotNsubgroups";
    case SUBFIND_IDS:
    case SNIP_SUBFIND_IDS:
      return "/Header/TotNids";
    case SNAP:
    case SNIP:
    case PARTDATA: 
    case SNIP_PARTDATA: 
      return "/Header/NumPart_Total";
    default:
      error("Type of file unknown");
    }
}

/**
 * \brief Returns the header field containg the number of elements in this file
 */
inline std::string elementFieldName(fileType type)
{
  switch(type)
    {
    case FOF_PARTICLES: 
    case SUBFIND_IDS:
    case SNIP_FOF_PARTICLES: 
    case SNIP_SUBFIND_IDS:
      return "/Header/Nids";
    case FOF: 
    case SUBFIND_GROUP:
    case SNIP_FOF:
    case SNIP_SUBFIND_GROUP:
      return "/Header/Ngroups";
    case SUBFIND: 
    case SNIP_SUBFIND: 
      return "/Header/Nsubgroups";
    case SNAP:
    case SNIP:
    case PARTDATA: 
    case SNIP_PARTDATA: 
      return "/Header/NumPart_ThisFile";
    default:
      error("Type of file unknown");
    }
}

/**
 * \brief Returns the name of this fileType
 */
inline std::string fileTypeName(fileType type)
{
  switch(type)
    {
    case FOF: return "FOF";
    case SNIP_FOF: return "FOF (snipshot)";
    case FOF_PARTICLES: return "FOF (particle info)";
    case SNIP_FOF_PARTICLES: return "FOF (particle info, snipshot)";
    case SUBFIND: return "Subfind";
    case SUBFIND_GROUP: return "Subfind (group tabs)";
    case SUBFIND_IDS: return "Subfind (IDs)";
    case SNIP_SUBFIND: return "Subfind (snipshot)";
    case SNIP_SUBFIND_GROUP: return "Subfind (group tabs, snipshot)";
    case SNIP_SUBFIND_IDS: return "Subfind (IDs, snipshot)";
    case SNIP: return "Snipshot";
    case SNAP: return "Snapshot";
    case PARTDATA: return "ParticleData";
    case SNIP_PARTDATA: return "ParticleData (snipshot)";
    default:
      error("Type of file unknown");
    }
}

inline int getGadgetType(fileType type, const std::string& arrayName)
{
  switch(type)
    {
    case FOF:
    case FOF_PARTICLES:
    case SUBFIND:
    case SUBFIND_GROUP:
    case SUBFIND_IDS:
    case SNIP_FOF:
    case SNIP_FOF_PARTICLES:
    case SNIP_SUBFIND:
    case SNIP_SUBFIND_GROUP:
    case SNIP_SUBFIND_IDS:
      return 0;
    case SNIP:
    case SNAP:
    case PARTDATA:
    case SNIP_PARTDATA:
      if(arrayName[0] != '/')
	error("Ill-formed array name. Must start with a /");
      return atoi(arrayName.substr(9,1).c_str());
    default:
      error("Type of file unknown");
    }
}


/**
 *\brief  Returns a vector of strings containing the different elements of a path of the form /aaa/bbb/ccc
 */
inline std::vector<std::string> tokenizePath(const std::string& path)
{
  std::vector<std::string> tokens;
  const char delimiter='/';
  std::string::size_type lastPos = path.find_first_not_of(delimiter, 0);
  std::string::size_type pos     = path.find_first_of(delimiter, lastPos);
  while (pos != std::string::npos || lastPos !=std::string::npos)
    {
      tokens.push_back(path.substr(lastPos, pos - lastPos));
      lastPos = path.find_first_not_of(delimiter, pos);
      pos = path.find_first_of(delimiter, lastPos);
    }
  return tokens;
}

/**
 * \brief Returns the value of the attribute 'name' in an open HDF5 group.
 */
template<typename T>
inline T readAttributeFromGroup(hid_t grp, const std::string& name)
{
  T data;

  const hid_t h_attr = H5Aopen(grp, name.c_str(), H5P_DEFAULT);
  if(h_attr < 0)
    error("Impossible to open attribute");

  const herr_t h_err = H5Aread(h_attr, hdf5Type<T>(), &data);
  if(h_err < 0)
    error("Impossible to read attribute");

  H5Aclose(h_attr);
  return data;
}


/**
 * \brief Returns the value of the attribute 'name' in an open HDF5 file.
 */
template<typename T>
T readAttribute(hid_t h_file, const std::string& name)
{
  const std::vector<std::string> tokens = tokenizePath(name);
  std::vector<hid_t> groups(tokens.size(), 0);

  groups[0] = h_file;
  for(size_t i(1); i<groups.size(); ++i)
    groups[i] = H5Gopen1(groups[i-1], (tokens[i-1]).c_str());
    
  const T attr = readAttributeFromGroup<T>(groups.back(), tokens.back().c_str());

  for(size_t i(1); i<groups.size(); ++i)
    H5Gclose(groups[i]);
  return attr;
}

/**
 * \brief Returns the value of the attribute 'name' in an HDF5 file given by its dir and tag.
 */
template<typename T>
inline T readAttribute(fileType type, const std::string& dir, const std::string& tag, const std::string& name)
{
 
  const std::string fname=fileName(type, dir, tag);
  const hid_t h_file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  const T attr = readAttribute<T>(h_file, name);
  H5Fclose(h_file);
  return attr;
}



/**
 * \brief Returns the type of the attribute 'name' in an open HDF5 group.
 */
inline hid_t readAttributeTypeFromGroup(hid_t grp, const std::string& name)
{
  const hid_t h_attr = H5Aopen(grp, name.c_str(), H5P_DEFAULT);
  if(h_attr < 0)
    error("Impossible to open attribute");

  const herr_t h_type = H5Aget_type(h_attr);
  if(h_type < 0)
    error("Impossible to read attribute type");

  H5Aclose(h_attr);
  return h_type;
}


/**
 * \brief Returns the type of the attribute 'name' in an open HDF5 file.
 */
inline hid_t readAttributeType(hid_t h_file, const std::string& name)
{
  const std::vector<std::string> tokens = tokenizePath(name);
  std::vector<hid_t> groups(tokens.size(), 0);

  groups[0] = h_file;
  for(size_t i(1); i<groups.size(); ++i)
    groups[i] = H5Gopen1(groups[i-1], (tokens[i-1]).c_str());
    
  const hid_t attr = readAttributeTypeFromGroup(groups.back(), tokens.back().c_str());

  for(size_t i(1); i<groups.size(); ++i)
    H5Gclose(groups[i]);
  return attr;
}


/**
 * \brief Returns the HDF5 type of a given attribute 'name' in an HDF5 file given its dir and tag.
 */
inline hid_t readAttributeType(fileType type, const std::string& dir, const std::string& tag, const std::string& name)
{
  const std::string fname=fileName(type, dir, tag);
  const hid_t h_file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  const hid_t h_type = readAttributeType(h_file, name);
  H5Fclose(h_file);
  return h_type;
}

/**
 * \brief Returns the size of the attribute 'name' in an open HDF5 group.
 */
inline hid_t readAttributeSizeFromGroup(hid_t grp, const std::string& name)
{
  const hid_t h_attr = H5Aopen(grp, name.c_str(), H5P_DEFAULT);
  if(h_attr < 0)
    error("Impossible to open attribute");

  const hid_t h_space = H5Aget_space(h_attr);
  if(h_space < 0)
    error("Impossible to open attribute space");

  const htri_t isSimple = H5Sis_simple(h_space);
  if(isSimple < 0)
    error("Impossible to extract information about the attribute's dataspace");
  else if(isSimple == 0)
    error("Attribute's dataspace is not simple");
    
  hsize_t dims[1];
  const int ndims = H5Sget_simple_extent_dims(h_space, dims, 0);
  if(ndims < 0)
    error("Attribute's dataspace is weird");
  else if(ndims == 0)
    dims[0] = 1;

  H5Sclose(h_space);
  H5Aclose(h_attr);
  return dims[0];
}


/**
 * \brief Returns the size of the attribute 'name' in an open HDF5 file.
 */
inline int readAttributeSize(hid_t h_file, const std::string& name)
{
  const std::vector<std::string> tokens = tokenizePath(name);
  std::vector<hid_t> groups(tokens.size(), 0);

  groups[0] = h_file;
  for(size_t i(1); i<groups.size(); ++i)
    groups[i] = H5Gopen1(groups[i-1], (tokens[i-1]).c_str());
    
  const int size = readAttributeSizeFromGroup(groups.back(), tokens.back().c_str());

  for(size_t i(1); i<groups.size(); ++i)
    H5Gclose(groups[i]);
  return size;
}

/**
 * \brief Returns the size of a given attribute 'name' in an HDF5 file given its dir and tag.
 */
inline int readAttributeSize(fileType type, const std::string& dir, const std::string& tag, const std::string& name)
{
  const std::string fname=fileName(type, dir, tag);
  const hid_t h_file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  const int size = readAttributeSize(h_file, name);
  H5Fclose(h_file);
  return size;
}


/**
 * \brief Returns the value of the attribute 'name' in an open HDF5 file.
 */
template<typename T>
T readUnits(hid_t h_file, const std::string& name)
{
  const std::vector<std::string> tokens = tokenizePath(name);
  std::vector<hid_t> groups(tokens.size(), 0);

  groups[0] = h_file;
  for(size_t i(1); i<groups.size()-1; ++i)
    groups[i] = H5Gopen1(groups[i-1], (tokens[i-1]).c_str());
    
  const hid_t h_data = H5Dopen1(groups[groups.size()-2], tokens[groups.size()-2].c_str());
  const T attr = readAttributeFromGroup<T>(h_data, tokens.back().c_str());
  H5Dclose(h_data);

  for(size_t i(1); i<groups.size()-1; ++i)
    H5Gclose(groups[i]);
  return attr;
}



/**
 * \brief Returns the type of the attribute 'name' in an open HDF5 group.
 */
inline hid_t readArrayTypeFromGroup(hid_t grp, const std::string& name)
{
  const hid_t h_attr = H5Dopen(grp, name.c_str(), H5P_DEFAULT);
  if(h_attr < 0)
    error("Impossible to open array");

  const herr_t h_type = H5Dget_type(h_attr);
  if(h_type < 0)
    error("Impossible to read array type");

  H5Dclose(h_attr);
  return h_type;
}


/**
 * \brief Returns the type of the attribute 'name' in an open HDF5 file.
 */
inline hid_t readArrayType(hid_t h_file, const std::string& name)
{
  const std::vector<std::string> tokens = tokenizePath(name);
  std::vector<hid_t> groups(tokens.size(), 0);

  // for(int i(0); i<tokens.size();++i)
  //   std::cout << tokens[i] << std::endl;

  groups[0] = h_file;
  for(size_t i(1); i<groups.size(); ++i)
    {
      // std::cout << "Opening group '" << tokens[i-1] << "'" <<std::endl;
      groups[i] = H5Gopen1(groups[i-1], (tokens[i-1]).c_str());
      // std::cout << "Group open" << std::endl;
    }
    
  const hid_t attr = readArrayTypeFromGroup(groups.back(), tokens.back().c_str());

  for(size_t i(1); i<groups.size(); ++i)
    H5Gclose(groups[i]);
  return attr;
}


/**
 * \brief Returns the HDF5 type of a given attribute 'name' in an HDF5 file given its dir and tag.
 */
inline hid_t readArrayType(fileType type, const std::string& dir, const std::string& tag, const std::string& name)
{
  const std::string fname=fileName(type, dir, tag);
  const hid_t h_file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  const hid_t h_type = readArrayType(h_file, name);
  H5Fclose(h_file);
  return h_type;
}

/**
 * \brief Returns the size of the array 'name' in an open HDF5 group.
 */
inline hid_t readArraySizeFromGroup(hid_t grp, const std::string& name)
{
  const hid_t h_attr = H5Dopen(grp, name.c_str(), H5P_DEFAULT);
  if(h_attr < 0)
    error("Impossible to open dataset");

  const hid_t h_space = H5Dget_space(h_attr);
  if(h_space < 0)
    error("Impossible to open data space");

  const htri_t isSimple = H5Sis_simple(h_space);
  if(isSimple < 0)
    error("Impossible to extract information about the array's dataspace");
  else if(isSimple == 0)
    error("array's dataspace is not simple");
    
  hsize_t dims[2];
  const int ndims = H5Sget_simple_extent_dims(h_space, dims, 0);
  if(ndims < 0)
    error("Array's dataspace is weird");

  H5Sclose(h_space);
  H5Dclose(h_attr);
  if(ndims == 1)
    return 1;
  else if(ndims == 2)
    return dims[1];
  else
    error("Array's dataspace is weird");
}


/**
 * \brief Returns the size of the array 'name' in an open HDF5 file.
 */
inline int readArraySize(hid_t h_file, const std::string& name)
{
  const std::vector<std::string> tokens = tokenizePath(name);
  std::vector<hid_t> groups(tokens.size(), 0);

  groups[0] = h_file;
  for(size_t i(1); i<groups.size(); ++i)
    groups[i] = H5Gopen1(groups[i-1], (tokens[i-1]).c_str());
    
  const int size = readArraySizeFromGroup(groups.back(), tokens.back().c_str());

  for(size_t i(1); i<groups.size(); ++i)
    H5Gclose(groups[i]);
  return size;
}

/**
 * \brief Returns the size of a given array 'name' in an HDF5 file given its dir and tag.
 */
inline int readArraySize(fileType type, const std::string& dir, const std::string& tag, const std::string& name)
{
  const std::string fname=fileName(type, dir, tag);
  const hid_t h_file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  const int size = readArraySize(h_file, name);
  H5Fclose(h_file);
  return size;
}



template<typename T, typename S> 
void readArray(hid_t file, const std::string& name, S* data, size_t offset=0)
{
  const hid_t h_data = H5Dopen1(file, name.c_str());
  if(h_data < 0)
    error("Impossible to open file");

  const hid_t h_space = H5Dget_space(h_data);
  if(h_space < 0)
    error("Impossible to open dataspace");

  hsize_t dims[2], dummy[2];
  int rank = H5Sget_simple_extent_dims(h_space, dims, dummy);

  if(rank == 2)
    offset *= dims[1];

  const hid_t h_err = H5Dread(h_data, hdf5Type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, reinterpret_cast<T*>(data) + offset);
  if(h_err < 0)
    error("Impossible to read dataset");

  H5Dclose(h_data);
}


/**
 * \brief Reads the content of the dataset 'name' into the vector 'data' in an open HDF5 file.
 */
template<typename T, typename S> 
inline void readArray(hid_t file, const std::string& name, std::vector<S>& data, size_t offset=0)
{
  readArray(file, name, &data[0], offset);
}

template <typename S>
struct thread_info
{
  int id;
  int numSnap;
  int numThreads;
  fileType type;
  int partType;
  std::string dir;
  std::string tag;
  std::string arrayName;
  S* array;
  size_t* sizes;
  size_t* offsets;
  bool oldSubfindFix;
};

template<typename S>
void* thread_readSize(void *thread_data)
{
  thread_info<S> t_data = *reinterpret_cast<thread_info<S>*>( thread_data );

  for(int i(t_data.id); i<t_data.numSnap; i+=t_data.numThreads)
    {
      //std::cerr << "Thread " << t_data.id << "/" << t_data.numThreads << " reading file " << i << "/" << t_data.numSnap  << std::endl;
      const hid_t file = H5Fopen(fileName(t_data.type, t_data.dir, t_data.tag, i).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if(file == 0)
	error("Impossible to open file");

      if(!t_data.oldSubfindFix)
	t_data.sizes[i] = readAttribute<Tuple<6, int> >(file, elementFieldName(t_data.type))[t_data.partType];
      else
	{
	  t_data.sizes[i] = readAttribute<Tuple<6, int> >(file, "/Header/NumPart_ThisFile")[2];
	}

      H5Fclose(file);
    }
  return 0;
}

template<typename T, typename S>
void* thread_readData(void *thread_data)
{
  thread_info<S> t_data = *reinterpret_cast<thread_info<S>*>( thread_data );

  for(int i(t_data.id); i<t_data.numSnap; i+=t_data.numThreads)
    {
      // std::cerr << "Thread " << t_data.id << "/" << t_data.numThreads << " reading file " << i << "/" << t_data.numSnap << " " << t_data.arrayName << " " <<  t_data.array << " " <<  t_data.offsets[i] << " " << t_data.type << std::endl;
      const hid_t file = H5Fopen(fileName(t_data.type, t_data.dir, t_data.tag, i).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if(file == 0)
	error("Impossible to open file");
      if(t_data.sizes[i] > 0)
	{
	  readArray<T>(file, t_data.arrayName, t_data.array, t_data.offsets[i]);
	}
      H5Fclose(file);

      // std::cerr << "Done" << std::endl;
    }

  return 0;
}

/**
 * \brief Reads the content of the dataset 'arrayName' into the vector 'array' for a given file (type, dir, tag)
 */
template<typename T, typename S>
void readArray(fileType type, const std::string& dir, const std::string& tag, const std::string& arrayName, std::vector<S>& array, bool physicalUnits=true, bool noH=true, int numThreads=16, bool verbose=true)
{
  error("Strange function call");

  /* Check that we have at least one thread */
  if(numThreads < 1)
    error("Invalid number of threads. Value must be >0.");

 /* Open file #0 to read global information */  
  const hid_t h_file = H5Fopen(fileName(type, dir, tag).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if(h_file < 0)
    error("Impossible to read file #0");

  /* Decide which of the 6 values for the number of elements to read if we are dealing with a sn[ai]pshot */
  const int partType = getGadgetType(type, arrayName);

  /* Decide which of the 6 values for the number of elements to read if we are dealing with a sn[ai]pshot */
  const int numSnap = readAttribute<int>(h_file, numFilesFieldName(type));
  const size_t totNumElem = readAttribute<Tuple<6, int> >(h_file, totElementFieldName(type))[partType];
  const double hubbleParam = readAttribute<double>(h_file, "/Header/HubbleParam");
  const double a = readAttribute<double>(h_file, "/Header/ExpansionFactor");
  H5Fclose(h_file);
 
  if(verbose)
    std::cout << "Extracting " << totNumElem << " values for '"<< arrayName << "' (PartType=" << partType << ") in " << numSnap << " " << fileTypeName(type) << " files with tag '" << tag << "' using " << numThreads << " thread(s)." << std::endl;
  
  /* Prepare array */
  array.resize(totNumElem);
  if(totNumElem == 0)
    return;

  /* Separate arrays for parallel reads */
  std::vector<size_t> sizes(numSnap, 0);
  std::vector<size_t> offsets(numSnap, 0);
  
  /* Prepare data for threads*/
  pthread_t* threads = new pthread_t[numThreads];
  thread_info<S>* thread_out = new thread_info<S>[numThreads];
  for(int i(0); i<numThreads; ++i)
    {
      thread_out[i].id = i;
      thread_out[i].numSnap = numSnap;
      thread_out[i].numThreads = numThreads;
      thread_out[i].type = type;
      thread_out[i].partType = partType;
      thread_out[i].dir = dir;
      thread_out[i].tag = tag;
      thread_out[i].arrayName = arrayName;
      thread_out[i].array = &array[0];
      thread_out[i].sizes = &sizes[0];
      thread_out[i].offsets = &offsets[0];
      thread_out[i].oldSubfindFix = false;
    }

  timeval begin, end;
  gettimeofday(&begin, NULL);

  if(verbose)
    std::cout << "Reading array lengths" << std::flush;

  /* Compute sizes for parallel reads using threads */
  for(int i(0); i<numThreads; ++i)
    {
      const int ret = pthread_create(&threads[i], 0, thread_readSize<S>, &thread_out[i]);
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


  /* Do a scan to convert from size to offset */
  for(int i(1); i<numSnap; ++i)
    {
      offsets[i] = offsets[i-1] + sizes[i-1];
      //      std::cout << i << " " << offsets[i] << " " << sizes[i] << std::endl;
    }


  /* Read data from files using threads */
  gettimeofday(&begin, NULL);
  if(verbose)
    std::cout << "Reading data" << std::flush;

  for(int i(0); i<numThreads; ++i)
    {
      const int ret = pthread_create(&threads[i], 0, thread_readData<T,S>, &thread_out[i]);
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
  //debug_variable(nonEmpty);

  /* Read factors for unit conversion */
  const hid_t h_file2 = H5Fopen(fileName(type, dir, tag, nonEmpty).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  const float a_exponant = readUnits<float>(h_file2, arrayName + "/aexp-scale-exponent");
  const float h_exponant = readUnits<float>(h_file2, arrayName + "/h-scale-exponent");
  H5Fclose(h_file2);

  /* Apply comoving to physical units conversion */
  if(physicalUnits)
    {
      if(a_exponant != 0.)
	{
	  if(verbose)
	    std::cout << "Converting to physical units. (Multiplication by a^" << a_exponant << ", a=" << a  << ")" << std::endl;
	  const double factor = std::pow(a, static_cast<double>(a_exponant));
	  for(size_t i(0); i<totNumElem; ++i)
	    array[i] *= factor;
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
	  for(size_t i(0); i<totNumElem; ++i)
	    array[i] *= factor;
	}
      else
	{
	  if(verbose)
	    std::cout << "Converting to h-free units. No conversion needed !" << std::endl;
	}
    }

  delete[] threads;
  delete[] thread_out;
}








template<typename T, typename S>
void readMultipleArrays(fileType type, const std::string& dir, const std::string& tag, const std::vector<std::string>& arrayName, std::vector<std::vector<S> >& array, bool physicalUnits=true, bool noH=true, int numThreads=16, bool verbose=true)
{
  /* Check that we have at least one thread */
  if(numThreads < 1)
    error("Invalid number of threads. Value must be >0.");

  /* Open file #0 to read global information */  
  const hid_t h_file = H5Fopen(fileName(type, dir, tag).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if(h_file < 0)
    error("Impossible to read file #0");

  const size_t numArrays = arrayName.size();
  if(numArrays != array.size())
    error("'arrayName' and 'array' have different sizes");

  /* Decide which of the 6 values for the number of elements to read if we are dealing with a sn[ai]pshot */
  const int partType = getGadgetType(type, arrayName[0]);

  /* Check that the partType is the same for all arrays */
  for(size_t i(1); i<arrayName.size(); ++i)
    if(partType != getGadgetType(type, arrayName[i]))
      error("Not all arrays have the same particle type");

  /* Read general information from file #0 */
  const int numSnap = readAttribute<int>(h_file, numFilesFieldName(type));
  const int totNumElem = readAttribute<Tuple<6, int> >(h_file, totElementFieldName(type))[partType];
  const double hubbleParam = readAttribute<double>(h_file, "/Header/HubbleParam");
  const double a = readAttribute<double>(h_file, "/Header/ExpansionFactor");
  H5Fclose(h_file);

  if(verbose)
    std::cout << "Extracting " << totNumElem << " values for many arrays (PartType=" << partType << ") in " << numSnap << " " << fileTypeName(type) << " files with tag '" << tag << "' using " << numThreads << " thread(s)." << std::endl;
  
  /* Prepare array */
  for(size_t i(0); i<numArrays; ++i)
    array[i].resize(totNumElem);
  if(totNumElem == 0)
    return;

  /* Separate arrays for parallel reads */
  std::vector<size_t> sizes(numSnap, 0);
  std::vector<size_t> offsets(numSnap, 0);
  
  /* Prepare data for threads*/
  pthread_t* threads = new pthread_t[numThreads];
  thread_info<S>* thread_out = new thread_info<S>[numThreads];
  for(int i(0); i<numThreads; ++i)
    {
      thread_out[i].id = i;
      thread_out[i].numSnap = numSnap;
      thread_out[i].numThreads = numThreads;
      thread_out[i].type = type;
      thread_out[i].partType = partType;
      thread_out[i].dir = dir;
      thread_out[i].tag = tag;
      thread_out[i].arrayName = arrayName[0];
      thread_out[i].array = &array[0][0];
      thread_out[i].sizes = &sizes[0];
      thread_out[i].offsets = &offsets[0];
      thread_out[i].oldSubfindFix = false;
    }

  timeval begin, end;
  gettimeofday(&begin, NULL);

  if(verbose)
    std::cout << "Reading array lengths" << std::flush;

  /* Compute sizes for parallel reads using threads */
  for(int i(0); i<numThreads; ++i)
    {
      const int ret = pthread_create(&threads[i], 0, thread_readSize<S>, &thread_out[i]);
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


  /* Do a scan to convert from size to offset */
  for(int i(1); i<numSnap; ++i)
    offsets[i] = offsets[i-1] + sizes[i-1];

  /* Loop over the arrays to read */
  for(size_t j(0); j<numArrays; ++j)
    {

      for(int i(0); i<numThreads; ++i)
	{
	  thread_out[i].arrayName = arrayName[j];
	  thread_out[i].array = &array[j][0];
	}

      /* Read data from files using threads */
      gettimeofday(&begin, NULL);
      if(verbose)
	std::cout << "Reading data (" << arrayName[j] << ") " << std::flush;

      for(int i(0); i<numThreads; ++i)
	{
	  const int ret = pthread_create(&threads[i], 0, thread_readData<T,S>, &thread_out[i]);
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

      /* Read factors for unit conversion */
      const hid_t h_file2 = H5Fopen(fileName(type, dir, tag, nonEmpty).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      const float a_exponant = readUnits<float>(h_file2, arrayName[j] + "/aexp-scale-exponent");
      const float h_exponant = readUnits<float>(h_file2, arrayName[j] + "/h-scale-exponent");
      H5Fclose(h_file2);

      /* Apply comoving to physical units conversion */
      if(physicalUnits && a_exponant != 0.)
	{
	  const double factor = std::pow(a, static_cast<double>(a_exponant));
	  for(int i(0); i<totNumElem; ++i)
	    array[j][i] *= factor;
	}
      
      /* Apply h units conversion */
      if(noH && h_exponant != 0.)
	{
	  const double factor = std::pow(hubbleParam, static_cast<double>(h_exponant));
	  for(int i(0); i<totNumElem; ++i)
	    array[j][i] *= factor;
	}
    }

  delete[] threads;
  delete[] thread_out;
}



#endif //READ_ARRAY_HPP
