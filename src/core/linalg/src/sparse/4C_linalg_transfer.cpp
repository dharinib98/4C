// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later



#include "4C_config.hpp"

#include "4C_linalg_transfer.hpp"

// Do not lint the file for identifier names, since the naming of the Wrapper functions follow the
// naming of the Core::LinAlg::Map

// NOLINTBEGIN(readability-identifier-naming)

FOUR_C_NAMESPACE_OPEN

// Standard constructor using source and target maps
Core::LinAlg::Import::Import::Import(
    const Epetra_BlockMap& targetMap, const Epetra_BlockMap& sourceMap)
    : import_(std::make_shared<Epetra_Import>(targetMap, sourceMap))
{
}

// Copy constructor
Core::LinAlg::Import::Import::Import(const Import& other)
    : import_(std::make_shared<Epetra_Import>(*other.import_))
{
}

// Construct from Epetra_Export (reverse communication)
Core::LinAlg::Import::Import::Import(const Epetra_Export& exporter)
    : import_(std::make_shared<Epetra_Import>(exporter))
{
}

// Expert constructor: pre-computed remote PIDs
Core::LinAlg::Import::Import::Import(const Epetra_BlockMap& targetMap,
    const Epetra_BlockMap& sourceMap, int numRemotePIDs, const int* remotePIDs)
    : import_(std::make_shared<Epetra_Import>(targetMap, sourceMap, numRemotePIDs, remotePIDs))
{
}

// Expert constructor: full manual setup
Core::LinAlg::Import::Import::Import(const Epetra_BlockMap& targetMap,
    const Epetra_BlockMap& sourceMap, int numRemotePIDs, const int* remotePIDs,
    const int& numExportIDs, const int* exportLIDs, const int* exportPIDs)
    : import_(std::make_shared<Epetra_Import>(
          targetMap, sourceMap, numRemotePIDs, remotePIDs, numExportIDs, exportLIDs, exportPIDs))
{
}

// Getter for raw Epetra_Import reference
const Epetra_Import& Core::LinAlg::Import::Import::get_epetra_import() const { return *import_; }

// Forwarded access to the underlying SourceMap
const Epetra_BlockMap& Core::LinAlg::Import::SourceMap() const { return import_->SourceMap(); }

// Forwarded access to the underlying TargetMap
const Epetra_BlockMap& Core::LinAlg::Import::TargetMap() const { return import_->TargetMap(); }

// NOLINTEND(readability-identifier-naming)

FOUR_C_NAMESPACE_CLOSE
