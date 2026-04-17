// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for output functions
//
// -----------------------------------------------------------------------------

#ifndef OUTPUT_UT_H
#define OUTPUT_UT_H

#include "unittest.h"
#include "IO/output.h"

// -----------------------------------------------------------------------------
// Description: Output UTs
// -----------------------------------------------------------------------------
void Test_output_writeModalSolution1DHeader();
void Test_output_writeModalSolution1DData();
void Test_output_writeSolution1D();
void Test_output_TimeSeriesWriterWriteSequence();

// -----------------------------------------------------------------------------
// Description: Output UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_output(TestRegistry& registry)
{
  registry.add("Test_output_writeModalSolution1DHeader",    Test_output_writeModalSolution1DHeader);
  registry.add("Test_output_writeModalSolution1DData",      Test_output_writeModalSolution1DData);
  registry.add("Test_output_writeSolution1D",               Test_output_writeSolution1D);
  registry.add("Test_output_TimeSeriesWriterWriteSequence", Test_output_TimeSeriesWriterWriteSequence);
}

#endif
