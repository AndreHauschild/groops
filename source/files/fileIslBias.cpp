/***********************************************/
/**
* @file fileIslBias.cpp
*
* @brief ISL biases.
*
* @author Andre Hauschild
* @date 2026-02-02
*
*/
/***********************************************/

#define DOCSTRING_FILEFORMAT_IslBias

#include "base/import.h"
#include "inputOutput/fileArchive.h"
#include "files/fileFormatRegister.h"
#include "fileIslBias.h"

GROOPS_REGISTER_FILEFORMAT(IslBias, FILE_ISLBIAS_TYPE)

/***********************************************/

Vector IslBias::compute(const std::vector<UInt> &terminals) const
{
  try
  {
    Vector b(terminals.size());
    UInt idx;
    for(UInt idTerm=0; idTerm<terminals.size(); idTerm++)
      if(isInList(terminals.at(idTerm), idx))
        b(idTerm) = biases.at(idx);
    return b;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

Bool IslBias::isInList(UInt terminal, UInt &idx) const
{
  idx = index(terminal);
  return (idx != NULLINDEX);
}

/***********************************************/

UInt IslBias::index(UInt terminal) const
{
  for(UInt i=0; i<terminals.size(); i++)
    if(terminal == terminals.at(i))
      return i;
  return NULLINDEX;
}

/***********************************************/
/***********************************************/

template<> void save(OutArchive &ar, const IslBias &x)
{
  try
  {
    ar<<nameValue("count", x.terminals.size());
    ar.comment(" terminal  bias [m]");
    ar.comment("===================================");
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

template<> void load(InArchive  &ar, IslBias &x)
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

void writeFileIslBias(const FileName &fileName, const IslBias &x)
{
  try
  {
    OutFileArchive file(fileName, FILE_ISLBIAS_TYPE, FILE_ISLBIAS_VERSION);
    file<<nameValue("islBias", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

void readFileIslBias(const FileName &fileName, IslBias &x)
{
  try
  {
    InFileArchive file(fileName, FILE_ISLBIAS_TYPE, FILE_ISLBIAS_VERSION);
    file>>nameValue("islBias", x);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
