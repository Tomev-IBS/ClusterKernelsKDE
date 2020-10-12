#ifndef CLUSTERKERNELSKDE_FUNCTION_H
#define CLUSTERKERNELSKDE_FUNCTION_H

#include <memory>

typedef std::vector<double> Point;

class RealValuedFunction{
  /*
   *  This class provides an interface for real valued functions of arbitrary input / output dimensions.
   *  It's especially important for this project, as Kernel is such a function.
   */
  public:
    virtual Point GetValue(const Point &pt) = 0;
};

typedef std::shared_ptr<RealValuedFunction> RealValuedFunctionPtr;

#endif //CLUSTERKERNELSKDE_FUNCTION_H
