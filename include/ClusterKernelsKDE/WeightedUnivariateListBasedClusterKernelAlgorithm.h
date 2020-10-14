#ifndef KERDEP_WEIGHTEDUNIVARIATELISTBASEDCLUSTERKERNELALGORITHM_H
#define KERDEP_WEIGHTEDUNIVARIATELISTBASEDCLUSTERKERNELALGORITHM_H

#include "UnivariateListBasedClusterKernelAlgorithm.h"

class WeightedUnivariateListBasedClusterKernelAlgorithm : public UnivariateListBasedClusterKernelAlgorithm {
  public:
    WeightedUnivariateListBasedClusterKernelAlgorithm(const int &m,
                                                      ClusterKernel*(*cluster_kernel_factory_method)(ClusterKernelStreamElement *stream_element));
    void PerformStep(ClusterKernelStreamElement *stream_element) override;
  protected:
    double weight_modifier_ = 0.2;
    void AddNewClusterKernel(ClusterKernelStreamElement *stream_element) override;
};

#endif //KERDEP_WEIGHTEDUNIVARIATELISTBASEDCLUSTERKERNELALGORITHM_H
