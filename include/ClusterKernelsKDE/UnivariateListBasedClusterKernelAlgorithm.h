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
    void UpdateMergeCostsListAfterAddingKernel(const unsigned int &new_kernel_position);
    void MergeClusterKernelsWithTheLowestMergeCost() override;
};

#endif //CLUSTERKERNELSKDE_UNIVARIATELISTBASEDCLUSTERKERNELALGORITHM_H
