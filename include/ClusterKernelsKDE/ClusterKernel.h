#ifndef CLUSTERKERNELSKDE_CLUSTERKERNEL_H
#define CLUSTERKERNELSKDE_CLUSTERKERNEL_H

#include <vector>
#include <memory>

#include "ClusterKernelStreamElement.h"

typedef std::vector<double> point;

class ClusterKernel{
  public:
    virtual point GetMean() = 0;
    virtual void Update(ClusterKernelStreamElement *stream_element) = 0;
    ClusterKernel* Merge(const ClusterKernel &other_cluster_kernel);
};

typedef std::shared_ptr<ClusterKernel> ClusterKernelPointer;

#endif //CLUSTERKERNELSKDE_CLUSTERKERNEL_H
