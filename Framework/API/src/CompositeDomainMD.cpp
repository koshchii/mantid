// Mantid Repository : https://github.com/mantidproject/mantid
//
// Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
//   NScD Oak Ridge National Laboratory, European Spallation Source,
//   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
// SPDX - License - Identifier: GPL - 3.0 +
//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidAPI/CompositeDomainMD.h"
#include "MantidAPI/FunctionDomainMD.h"
#include "MantidAPI/IMDIterator.h"

#include <stdexcept>

namespace Mantid::API {

/**
 * Create a composite domain from a IMDWorkspace.
 * @param ws :: Pointer to a workspace.
 * @param maxDomainSize :: The maximum size each domain can have.
 */
CompositeDomainMD::CompositeDomainMD(const IMDWorkspace_const_sptr &ws, size_t maxDomainSize)
    : m_iterator(ws->createIterator()) {
  m_totalSize = m_iterator->getDataSize();

  size_t maxDomainSizeDiv = maxDomainSize + 1;
  size_t nParts = m_totalSize / maxDomainSizeDiv;

  m_domains.resize(nParts);
  for (size_t i = 0; i < nParts - 1; ++i) {
    size_t start = i * maxDomainSize;
    m_domains[i] = std::make_unique<FunctionDomainMD>(ws, start, maxDomainSize);
  }
  size_t start = (nParts - 1) * maxDomainSize;
  m_domains.back() = std::make_unique<FunctionDomainMD>(ws, start, m_totalSize - start);
}

/**
 * Destructor.
 */
CompositeDomainMD::~CompositeDomainMD() = default;

/// Return i-th domain reset to its start.
const FunctionDomain &CompositeDomainMD::getDomain(size_t i) const {
  if (i >= m_domains.size()) {
    throw std::out_of_range("Domain index out of range");
  }
  m_domains[i]->reset();

  return *m_domains[i];
}

} // namespace Mantid::API
