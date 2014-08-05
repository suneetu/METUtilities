#include "GaudiKernel/DeclareFactoryEntries.h"

// Top level tool
#include "METUtilities/METRebuilder.h"
// Algs
#include "../METUtilAlg.h"

using namespace met;

DECLARE_TOOL_FACTORY(METRebuilder)
//
DECLARE_ALGORITHM_FACTORY(METUtilAlg)

DECLARE_FACTORY_ENTRIES(METReconstruction) {
  DECLARE_TOOL(METRebuilder)
  //
  DECLARE_ALGORITHM(METUtilAlg)
}
