// -----------------------------------------------------------------------------
//              DGmini, a minimal 1D discontinuous Galerkin solver              
// -----------------------------------------------------------------------------
//
// Description: Unittest utilities for testing DGmini building blocks
//
// -----------------------------------------------------------------------------

#ifndef UNITTEST_H
#define UNITTEST_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Description: UT name and function pointer
// -----------------------------------------------------------------------------
struct Test
{
  std::string name;
  void (*func)();
};

// -----------------------------------------------------------------------------
// Description: UT registry
// -----------------------------------------------------------------------------
class TestRegistry
{
  public:
    // -------------------------------------------------------------------------
    // Construction
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Access
    // -------------------------------------------------------------------------
    const std::vector<Test>& tests() const
    {
      return tests_;
    }

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------
    void add(const std::string& name, void (*func)())
    {
      tests_.push_back({name, func});
    }

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    
  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    std::vector<Test> tests_;
};

// -----------------------------------------------------------------------------
// Description: Boolean check
// -----------------------------------------------------------------------------
inline void Check(bool condition, const std::string& message)
{
  if (!condition)
  {
    throw std::runtime_error
    (
      message
  );
  }
}

// -----------------------------------------------------------------------------
// Description: Floating-point comparison
// -----------------------------------------------------------------------------
inline void CheckEqual(double actual,
                        double expected,
                        double tol,
                        const std::string& message)
{
  if (std::abs(actual - expected) > tol)
  {
    throw std::runtime_error
    (
      message + " | expected = " + std::to_string(expected) 
              +  ", actual = "   + std::to_string(actual)
    );
  }
}

#endif