/***********************************************/
/**
* @file fileIslSignalBias.cpp
*
* @brief ISL biases.
*
* @author Andre Hauschild
* @date 2026-02-02
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_IslSignalBias

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "files/fileIslSignalBias.h"

GROOPS_REGISTER_FILEFORMAT(IslSignalBias, FILE_ISLSIGNALBIAS_TYPE)

/***********************************************/

Double IslSignalBias::compute(const UInt terminal) const
{
  try
  {
    for(UInt idx=0; idx<this->terminals.size(); idx++)
      if(this->terminals.at(idx) == terminal)
        return biases.at(idx);
    return 0.0;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const IslSignalBias &x)
{
  try
  {
    ar<<nameValue("count", x.terminals.size());
    ar.comment("terminal bias [m]");
    ar.comment("===============================");
    for(UInt i=0; i<x.terminals.size(); i++)
    {
      ar.endLine();
      ar<<nameValue("terminal", x.terminals.at(i));
      ar<<nameValue("bias",     x.biases.at(i));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

template<> void load(InArchive  &ar, IslSignalBias &x)
{
  try
  {
    UInt count;
    ar>>nameValue("count", count);
    x.terminals.resize(count);
    x.biases.resize(count);
    for(UInt i=0; i<count; i++)
    {
      ar>>nameValue("terminal", x.terminals.at(i));
      ar>>nameValue("bias",     x.biases.at(i));
    }
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
/***********************************************/

void writeFileIslSignalBias(const FileName &fileName, const IslSignalBias &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_ISLSIGNALBIAS_TYPE, FILE_ISLSIGNALBIAS_VERSION);
    file<<nameValue("signalBias", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileIslSignalBias(const FileName &fileName, IslSignalBias &x)
{
  try
  {
    InFileArchive file(fileName, FILE_ISLSIGNALBIAS_TYPE, FILE_ISLSIGNALBIAS_VERSION);
    file>>nameValue("signalBias", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
