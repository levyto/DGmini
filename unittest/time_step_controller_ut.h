// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittests for TimeStepController class
//
// -----------------------------------------------------------------------------

#ifndef TIME_STEP_CONTROLLER_UT_H
#define TIME_STEP_CONTROLLER_UT_H

#include "unittest.h"
#include "Temporal/time_step_controller.h"

// -----------------------------------------------------------------------------
// Description: TimeStepController UTs
// -----------------------------------------------------------------------------
void Test_TimeStepController_FixedTimeStepReturnsGivenValue();
void Test_TimeStepController_CFLTimeStepLinearAdvectionConstantSolution();

// -----------------------------------------------------------------------------
// Description: TimeStepController UTs registry
// -----------------------------------------------------------------------------
inline void Register_Test_TimeStepController(TestRegistry& registry)
{
  registry.add("Test_TimeStepController_FixedTimeStepReturnsGivenValue",
                Test_TimeStepController_FixedTimeStepReturnsGivenValue);
  registry.add("Test_TimeStepController_CFLTimeStepLinearAdvectionConstantSolution",
                Test_TimeStepController_CFLTimeStepLinearAdvectionConstantSolution);
}
 
#endif