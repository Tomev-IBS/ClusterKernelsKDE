#include "WeightedUnivariateListBasedClusterKernelAlgorithm.h"

WeightedUnivariateListBasedClusterKernelAlgorithm::WeightedUnivariateListBasedClusterKernelAlgorithm(const int &m,
                                                                                                     ClusterKernel *(*cluster_kernel_factory_method)(
                                                                                                         ClusterKernelStreamElement *)) : UnivariateListBasedClusterKernelAlgorithm(m, cluster_kernel_factory_method)
{}

void WeightedUnivariateListBasedClusterKernelAlgorithm::PerformStep(ClusterKernelStreamElement *stream_element) {
  for(auto ck : cluster_kernels_){
    ck->RescaleWeight(1.0 - weight_modifier_);
  }
  ClusterKernelsAlgorithm::PerformStep(stream_element);
}

void WeightedUnivariateListBasedClusterKernelAlgorithm::AddNewClusterKernel(
    ClusterKernelStreamElement *stream_element) {
  auto new_cluster_kernel = ClusterKernelPointer(cluster_kernel_factory_method_(stream_element));
  new_cluster_kernel->RescaleWeight(weight_modifier_);
  // This is for univariate case. The idea is to place new cluster kernel in between
  // means closest to it.
  double new_cluster_kernel_value = new_cluster_kernel->GetMean()[0];
  unsigned int new_kernel_position = 0;
  for(auto value : cluster_kernels_){
    if(value->GetMean()[0] > new_cluster_kernel_value) {
      break;
    }
    ++new_kernel_position;
  }
  new_cluster_kernel->SetBandwidth(bandwidth_);
  cluster_kernels_.insert(cluster_kernels_.begin() + new_kernel_position, new_cluster_kernel);
  FillDomainForClusterKernelDistanceCalculation();
  UpdateMergeCostsListAfterAddingKernel(new_kernel_position);
}
