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
#include <iostream>

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
    int checks()   const { return checks_;   }
    int failures() const { return failures_; }

    // -------------------------------------------------------------------------
    // Modification
    // -------------------------------------------------------------------------
    void increment_checks()   { checks_++; }
    void increment_failures() { failures_++; }
    void add(const std::string& name, void (*func)())
    {
      tests_.push_back({name, func});
    }

    // -------------------------------------------------------------------------
    // Operations
    // -------------------------------------------------------------------------
    void run_all()
    {
      for (const auto& [name, test] : tests_)
      {
        const int failures_before = failures_;

        try
        {
          test();
        }
        catch (const std::exception& e)
        {
          // Just in case a test throws an exception, 
          // we catch it to report the failure and continue with other tests
          increment_failures();
          std::cout << "[EXCEPTION] " << name << " : " << e.what() << '\n';
        }

        if (failures_ == failures_before)
        {
          std::cout << "[PASS] " << name << '\n';
        }
        else
        {
          std::cout << "[FAIL] " << name << '\n';
        }
      }
    }
    
  private:
    // -------------------------------------------------------------------------
    // Data
    // -------------------------------------------------------------------------
    struct Test // UT name and function pointer
    {
      std::string name;
      void (*func)();
    };
    std::vector<Test> tests_;
    int checks_ = 0;
    int failures_ = 0;
};

// -----------------------------------------------------------------------------
// Description: Create a global instance of TestRegistry
// -----------------------------------------------------------------------------
inline TestRegistry& GetTestRegistry()
{
  static TestRegistry instance;
  return instance;
}

// -----------------------------------------------------------------------------
// Description: Boolean check
// -----------------------------------------------------------------------------
inline void Check(bool condition, const std::string& message)
{
  auto& reg = GetTestRegistry();
  reg.increment_checks();

  if (!condition)
  {
    reg.increment_failures();

    std::cout << "[ERROR] " << message << '\n';
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
  auto& reg = GetTestRegistry();
  reg.increment_checks();
  
  if (std::abs(actual - expected) > tol)
  {
    reg.increment_failures();

    std::cout << "[ERROR] "       << message
              << " | expected = " << expected
              << ", actual = "    << actual   << '\n';
  }
}

#endif