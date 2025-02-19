/***********************************************/
/**
* @file archiveJson.h
*
* @brief Read/write archive files in JSON format.
*
* @author Torsten Mayer-Guerr
* @date 2025-01-22
*
*/
/***********************************************/

#ifndef __GROOPS_ARCHIVEJSON__
#define __GROOPS_ARCHIVEJSON__

#include "archive.h"
#include "parser/xml.h"

/** @addtogroup archiveGroup */
/// @{

/***** CLASS ***********************************/

/** @brief Write archive files in JSON format. */
class OutArchiveJson : public OutArchive
{
public:
  OutArchiveJson(std::ostream &_stream, const std::string &type, UInt version);
 ~OutArchiveJson();

  OutArchiveJson() = delete;
  OutArchiveJson(const OutArchiveJson &) = delete;
  OutArchiveJson &operator=(const OutArchiveJson &) = delete;

  ArchiveType archiveType() const override {return JSON;}

protected:
  void startTag(const std::string &name) override;
  void endTag  (const std::string &name) override;

  void save(const std::string &x) override;
  void save(const Int         &x) override;
  void save(const UInt        &x) override;
  void save(const Double      &x) override;
  void save(const Bool        &x) override;
  void save(const Time        &x) override;
  void save(const Angle       &x) override;
  void save(const Vector      &x) override;
  void save(const const_MatrixSlice  &x) override;
  void save(const SphericalHarmonics &x) override;
  void save(const Doodson     &x) override;
  void save(const GnssType    &x) override;

private:
  std::ostream           &stream;
  std::stack<XmlNodePtr> stack;
};

/***** CLASS ***********************************/

/** @brief Read archive files in JSON format. */
class InArchiveJson : public InArchive
{
  std::stack<XmlNodePtr> stack;
  std::string typeStr;
  UInt       _version;

public:
  InArchiveJson(std::istream &stream);
 ~InArchiveJson();

  ArchiveType archiveType() const override {return JSON;}
  std::string type()        const override {return typeStr;}
  UInt        version()     const override {return _version;}

protected:
  void startTag(const std::string &name) override;
  void endTag  (const std::string &name) override;

  void load(std::string &x) override;
  void load(Int      &x) override;
  void load(UInt     &x) override;
  void load(Double   &x) override;
  void load(Bool     &x) override;
  void load(Time     &x) override;
  void load(Angle    &x) override;
  void load(Vector   &x) override;
  void load(Matrix   &x) override;
  void load(SphericalHarmonics &x) override;
  void load(Doodson  &x) override;
  void load(GnssType &x) override;
};

/***********************************************/

/// @}

#endif
