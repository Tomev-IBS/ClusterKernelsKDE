#ifndef CLUSTERKERNELSKDE_CLUSTERKERNELSALGORITHM_H
#define CLUSTERKERNELSKDE_CLUSTERKERNELSALGORITHM_H

#include <vector>

#include "ClusterKernelStreamElement.h"
#include "ClusterKernel.h"

class ClusterKernelsAlgorithm {

  public:
    ClusterKernelsAlgorithm(const int &m, ClusterKernel*(*cluster_kernel_factory_method)(ClusterKernelStreamElement *stream_element));
    void PerformStep(ClusterKernelStreamElement *stream_element);
  private:
    int maximal_number_of_cluster_kernels_ = 0;
    std::vector<ClusterKernelPointer> cluster_kernels_;
    ClusterKernel* (*cluster_kernel_factory_method_)(ClusterKernelStreamElement *stream_element);

    int FindIndexOfClusterKernelWithSameMeanAsStreamElement(ClusterKernelStreamElement *stream_element) const;
    double CalculateDistanceBetweenPoints(point point1, point point2) const;
    void AddNewClusterKernel(ClusterKernelStreamElement *stream_element);
    void MergeClusterKernelsWithTheLowestMergeCost();
};

#endif //CLUSTERKERNELSKDE_CLUSTERKERNELSALGORITHM_H
