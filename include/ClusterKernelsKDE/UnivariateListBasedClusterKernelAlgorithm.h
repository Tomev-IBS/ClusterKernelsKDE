#ifndef CLUSTERKERNELSKDE_UNIVARIATELISTBASEDCLUSTERKERNELALGORITHM_H
#define CLUSTERKERNELSKDE_UNIVARIATELISTBASEDCLUSTERKERNELALGORITHM_H

#include "ClusterKernelsAlgorithm.h"

class UnivariateListBasedClusterKernelAlgorithm : public ClusterKernelsAlgorithm {
  public:
    UnivariateListBasedClusterKernelAlgorithm(const int &m,
                                              ClusterKernel*(*cluster_kernel_factory_method)(ClusterKernelStreamElement *stream_element));
  protected:
    std::vector<double> merge_cost_with_next_cluster_kernel_;

    void AddNewClusterKernel(ClusterKernelStreamElement *stream_element) override;
    void UpdateMergeCostsList();
    double CalculateDistanceBetweenClusterKernelAndTheirMerge(const int &first_ck_index, const int &second_ck_index) override;
    void MergeClusterKernelsWithTheLowestMergeCost() override;
    void FillDomainForClusterKernelDistanceCalculation() override;
    double FindMinimalValueOnDimension(const int &dimension=0);
    double FindMaximalValueOnDimension(const int &dimension=0);
};

#endif //CLUSTERKERNELSKDE_UNIVARIATELISTBASEDCLUSTERKERNELALGORITHM_H
