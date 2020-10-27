#ifndef KERDEP_WEIGHTEDUNIVARIATELISTBASEDCLUSTERKERNELALGORITHM_H
#define KERDEP_WEIGHTEDUNIVARIATELISTBASEDCLUSTERKERNELALGORITHM_H

#include "UnivariateListBasedClusterKernelAlgorithm.h"

class WeightedUnivariateListBasedClusterKernelAlgorithm : public UnivariateListBasedClusterKernelAlgorithm {
  public:
    WeightedUnivariateListBasedClusterKernelAlgorithm(const int &m,
                                                      ClusterKernel*(*cluster_kernel_factory_method)(ClusterKernelStreamElement *stream_element));
    void PerformStep(ClusterKernelStreamElement *stream_element) override;
    Point GetValue(const Point &pt) override;
  protected:
    double weight_modifier_ = 0.01;
    double weights_sum_ = 0;
    void AddNewClusterKernel(ClusterKernelStreamElement *stream_element) override;
    double CalculateDistanceBetweenClusterKernelAndTheirMerge(const int &first_ck_index, const int &second_ck_index) override;
    void UpdateVariationEstimator(ClusterKernelStreamElement *stream_element) override;
};

#endif //KERDEP_WEIGHTEDUNIVARIATELISTBASEDCLUSTERKERNELALGORITHM_H
