// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_TRANSFER_HPP
#define FOUR_C_LINALG_TRANSFER_HPP

#include "4C_config.hpp"

#include <Epetra_Import.h>
#include <Epetra_Map.h>

#include <memory>

FOUR_C_NAMESPACE_OPEN

// Do not lint the file for identifier names, since the naming of the Wrapper functions follow the
// naming of the Epetra_Map

// NOLINTBEGIN(readability-identifier-naming)
namespace Core::LinAlg
{
  class Import
  {
   public:
    Import(const Epetra_BlockMap& targetMap, const Epetra_BlockMap& sourceMap);
    Import(const Import& other);
    Import(const Epetra_Export& exporter);
    Import(const Epetra_BlockMap& targetMap, const Epetra_BlockMap& sourceMap, int numRemotePIDs,
        const int* remotePIDs);
    Import(const Epetra_BlockMap& targetMap, const Epetra_BlockMap& sourceMap, int numRemotePIDs,
        const int* remotePIDs, const int& numExportIDs, const int* exportLIDs,
        const int* exportPIDs);

    const Epetra_Import& get_epetra_import() const;
    const Epetra_BlockMap& SourceMap() const;
    const Epetra_BlockMap& TargetMap() const;


   private:
    std::shared_ptr<Epetra_Import> import_;
  };
}  // namespace Core::LinAlg

// NOLINTEND(readability-identifier-naming)

FOUR_C_NAMESPACE_CLOSE
#endif
