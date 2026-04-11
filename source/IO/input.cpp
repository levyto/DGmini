// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definitions of inputs and readers
//
// -----------------------------------------------------------------------------

#include <stdexcept>

#include "IO/input.h"
#include "PDE/linear_advection1d.h"
#include "PDE/burgers1d.h"

// -----------------------------------------------------------------------------
// Description: Create a PDE instance based on its name
// -----------------------------------------------------------------------------
std::unique_ptr<PDE> createPDE(const std::string& name)
{
  if (name == "linear_advection1d")
    return std::make_unique<LinearAdvection1D>(1.0);

  if (name == "burgers")
    return std::make_unique<Burgers1D>();

  throw std::runtime_error("Unknown PDE model: " + name);
}
