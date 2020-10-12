#ifndef CLUSTERKERNELSKDE_CLUSTERKERNEL_H
#define CLUSTERKERNELSKDE_CLUSTERKERNEL_H

#include <vector>
#include <memory>

#include "ClusterKernelStreamElement.h"
#include "RealValuedFunction.h"

class ClusterKernel : public RealValuedFunction{
  public:
    virtual Point GetMean() = 0;
    virtual void Update(ClusterKernelStreamElement *stream_element) = 0;
    virtual ClusterKernel* Merge(ClusterKernel *other_cluster_kernel) = 0;
    Point GetValue(const Point &pt) override = 0;
    virtual double GetWeight() = 0;
    virtual unsigned int GetCardinality() = 0;
};

typedef std::shared_ptr<ClusterKernel> ClusterKernelPointer;

#endif //CLUSTERKERNELSKDE_CLUSTERKERNEL_H
