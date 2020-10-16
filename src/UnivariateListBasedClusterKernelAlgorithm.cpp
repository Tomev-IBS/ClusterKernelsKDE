#include "UnivariateListBasedClusterKernelAlgorithm.h"

UnivariateListBasedClusterKernelAlgorithm::UnivariateListBasedClusterKernelAlgorithm(const int &m,
                                                                                     ClusterKernel *(*cluster_kernel_factory_method)(
                                                                                         ClusterKernelStreamElement *)) : ClusterKernelsAlgorithm(m, cluster_kernel_factory_method)
{}

/** Adds new cluster kernel. This method is based on List-Based approach to cluster kernels algorithm for
 * univariate data. New element x should be inserted in such way that it's mean m_x satisfies
 *
 * m_{x-1} <= m_{x} <= m_{x+1}
 *
 * This will enable us to use faster finding of cluster kernels to merge.
 *
 * @brief Adds new cluster kernel.
 * @param stream_element
 */
void UnivariateListBasedClusterKernelAlgorithm::AddNewClusterKernel(ClusterKernelStreamElement *stream_element) {
  auto new_cluster_kernel = ClusterKernelPointer(cluster_kernel_factory_method_(stream_element));
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
  UpdateMergeCostsListAfterAddingKernel(new_kernel_position);
}

/** Updates merge costs list after adding new cluster kernel. Worst case scenario it requires calculation
 *  of two merge costs (for added element and previous one). Lowest cost scenario requires no calculations.
 *
 * Note: merge_cost_with_next_cluster_kernel_ has an property, namely it has the length of cluster_kernels_ - 1
 * at all times,
 *
 * @brief Updates merge costs list after adding new cluster kernel.
 * @param new_kernel_position - Position of new kernel.
 */
void UnivariateListBasedClusterKernelAlgorithm::UpdateMergeCostsListAfterAddingKernel(
    const unsigned int &new_kernel_position) {
  FillDomainForClusterKernelDistanceCalculation();
  merge_cost_with_next_cluster_kernel_.clear();
  for(auto i = 0; i < cluster_kernels_.size() - 1; ++i){
    merge_cost_with_next_cluster_kernel_.push_back(
      CalculateDistanceBetweenClusterKernelAndTheirMerge(i, i + 1);
    );
  }
}

/** Merges two cluster with the lowest merge cost. It takes advantage of the list-based
 * implementation of the algorithm to find the cluster kernels to merge faster. It also updates
 * the merge costs list at the end.
 *
 * @brief Merges two cluster with the lowest merge cost.
 */
void UnivariateListBasedClusterKernelAlgorithm::MergeClusterKernelsWithTheLowestMergeCost() {
  // I'll use abbreviation ck for cluster kernels.
  unsigned int index_of_ck_with_lowest_merge_cost_with_next_ck = 0;
  double minimal_distance = merge_cost_with_next_cluster_kernel_[index_of_ck_with_lowest_merge_cost_with_next_ck];

  for(auto current_ck_index = 1; current_ck_index < cluster_kernels_.size() - 1; ++current_ck_index){
    if(merge_cost_with_next_cluster_kernel_[current_ck_index] < minimal_distance){
      minimal_distance = merge_cost_with_next_cluster_kernel_[current_ck_index];
      index_of_ck_with_lowest_merge_cost_with_next_ck = current_ck_index;
    }
  }

  MergeClusterKernels(index_of_ck_with_lowest_merge_cost_with_next_ck,
                      index_of_ck_with_lowest_merge_cost_with_next_ck + 1);
}