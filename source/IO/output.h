// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Definitions of outputs
//
// -----------------------------------------------------------------------------

#ifndef OUTPUT_H
#define OUTPUT_H

#include <string>
#include <vector>

#include "Algebra/Vec.h"
#include "Mesh/mesh1d.h"

void writeModalSolution1D(const std::string& filename,
                          const Mesh1D& mesh,
                          const std::vector<Vec>& coefficients);

void writeSolution1D(const std::string& filename,
                     const Vec& x, 
                     const Vec& u);

#endif

