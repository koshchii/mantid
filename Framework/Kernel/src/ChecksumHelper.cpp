// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//     NScD Oak Ridge National Laboratory, European Spallation Source
//     & Institut Laue - Langevin
// SPDX - License - Identifier: GPL - 3.0 +
#include "MantidKernel/ChecksumHelper.h"

#include <Poco/DigestStream.h>
#include <Poco/MD5Engine.h>
#include <Poco/SHA1Engine.h>
#include <boost/regex.hpp>

#include <fstream>
#include <sstream>

namespace {

// Unix line ending character
constexpr auto EOL_LF = "\n";
// Windows line ending sequence
constexpr auto EOL_CRLF = "\r\n";

/**
 * Create sha1 out of data and an optional header
 * @param data Contents as a string
 * @param header An optional string to prepend to the data
 */
std::string createSHA1(const std::string &data,
                       const std::string &header = "") {
  using Poco::DigestEngine;
  using Poco::DigestOutputStream;
  using Poco::SHA1Engine;

  SHA1Engine sha1;
  DigestOutputStream outstr(sha1);
  outstr << header << data;
  outstr.flush(); // to pass everything to the digest engine
  return DigestEngine::digestToHex(sha1.digest());
}

/**
 * Create sha1 out of data and an optional header
 * @param data Contents as a string
 * @param header An optional string to prepend to the data
 */
std::string createMD5(const std::string &data, const std::string &header = "") {
  using Poco::DigestEngine;
  using Poco::DigestOutputStream;
  using Poco::MD5Engine;

  MD5Engine sha1;
  DigestOutputStream outstr(sha1);
  outstr << header << data;
  outstr.flush(); // to pass everything to the digest engine
  return DigestEngine::digestToHex(sha1.digest());
}
} // namespace

namespace Mantid::Kernel::ChecksumHelper {

/** Creates a md5 checksum from a string
 * @param input The string to checksum
 * @returns a checksum string
 **/
std::string md5FromString(const std::string &input) { return createMD5(input); }

/** Creates a SHA-1 checksum from a string
 * @param input The string to checksum
 * @returns a checksum string
 **/
std::string sha1FromString(const std::string &input) {
  return createSHA1(input);
}

/** Creates a SHA-1 checksum from a file
 * @param filepath The path to the file
 * @param unixEOL If true convert all lineendings to Unix-style \n
 * @returns a checksum string
 **/
std::string sha1FromFile(const std::string &filepath,
                         const bool unixEOL = false) {
  if (filepath.empty())
    return "";
  return createSHA1(loadFile(filepath, unixEOL));
}

/** Creates a git checksum from a file (these match the git hash-object
 *command).
 * This works by reading in the file, converting all line endings into linux
 *style endings,
 * then the following is prepended to the file contents "blob
 *<content_length>\0",
 * the result is then ran through a SHA-1 checksum.
 * @param filepath The path to the file
 * @returns a checksum string
 **/
std::string gitSha1FromFile(const std::string &filepath) {
  if (filepath.empty())
    return "";
  const bool unixEOL(true);
  std::string contents = ChecksumHelper::loadFile(filepath, unixEOL);
  std::stringstream header;
  header << "blob " << contents.size() << '\0';
  return createSHA1(contents, header.str());
}

/**
 * Load contents of file into a string. The line endings are preserved
 * @param filepath Full path to the file to be opened
 * @param unixEOL If true convert all lineendings to Unix-style \n
 */
std::string loadFile(const std::string &filepath, const bool unixEOL = false) {

  std::ifstream filein(filepath.c_str(), std::ios::in | std::ios::binary);
  if (!filein)
    return "";

  std::string contents;
  filein.seekg(0, std::ios::end);
  contents.resize(filein.tellg());
  filein.seekg(0, std::ios::beg);
  filein.read(&contents[0], contents.size());
  filein.close();

  if (unixEOL) {
    // convert CRLF to LF
    // Old-style mac line endings are left alone by git:
    // https://github.com/git/git/blob/bfdd66e72fffd18235757bedbf355fd4176d6019/convert.c#L1425
    static boost::regex eol(EOL_CRLF);
    contents = boost::regex_replace(contents, eol, EOL_LF);
  }
  return contents;
}

} // namespace Mantid::Kernel::ChecksumHelper
