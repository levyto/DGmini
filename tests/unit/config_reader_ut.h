// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for ConfigReader class
//
// -----------------------------------------------------------------------------

#ifndef CONFIG_READER_UT_H
#define CONFIG_READER_UT_H

#include "unittest.h"
#include "IO/config_reader.h"

// -----------------------------------------------------------------------------
// Description: ConfigReader UTs
// -----------------------------------------------------------------------------
void Test_ConfigReader_incrementalRequiredKeys();
void Test_ConfigReader_invalidEnumeratedValues();
void Test_ConfigReader_invalidNumericValues();
void Test_ConfigReader_optionalDefaults();

// -----------------------------------------------------------------------------
// Description: ConfigReader UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_ConfigReader(TestRegistry& registry)
{
  registry.add("Test_ConfigReader_incrementalRequiredKeys", Test_ConfigReader_incrementalRequiredKeys);
  registry.add("Test_ConfigReader_invalidEnumeratedValues", Test_ConfigReader_invalidEnumeratedValues);
  registry.add("Test_ConfigReader_invalidNumericValues",    Test_ConfigReader_invalidNumericValues   );
  registry.add("Test_ConfigReader_optionalDefaults",        Test_ConfigReader_optionalDefaults       );
}

#endif
