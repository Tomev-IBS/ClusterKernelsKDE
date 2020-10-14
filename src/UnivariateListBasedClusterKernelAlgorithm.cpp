#include "UnivariateListBasedClusterKernelAlgorithm.h"

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

  cluster_kernels_.insert(cluster_kernels_.begin() + new_kernel_position, new_cluster_kernel);
  UpdateMergeCostsListAfterAddingKernel(new_kernel_position);
}

/** Updates merge costs list after adding new cluster kernel. Worst case scenario it requires calculation
 *  of two merge costs (for added element and previous one). Lowest cost scenario requires no calculations.
 *
 * Note: merge_cost_with_next_cluster_kernel_ has two properties:
 * - has the same length same as cluster_kernels_ at all times,
 * - it's last element should always be -1;
 *
 * @brief Updates merge costs list after adding new cluster kernel.
 * @param new_kernel_position - Position of new kernel.
 */
void UnivariateListBasedClusterKernelAlgorithm::UpdateMergeCostsListAfterAddingKernel(
    const unsigned int &new_kernel_position) {
  // First thing I want to do is to insert a value. If there's no next cluster kernel
  // in the list I'll insert -1 (which is clearly not a distance measure value).
  auto new_kernel_iterator = merge_cost_with_next_cluster_kernel_.begin() + new_kernel_position;
  if(new_kernel_position + 1 == cluster_kernels_.size()){
    merge_cost_with_next_cluster_kernel_.insert(new_kernel_iterator, -1);
  } else {
    merge_cost_with_next_cluster_kernel_.insert(new_kernel_iterator, CalculateDistanceBetweenClusterKernelAndTheirMerge(new_kernel_position, new_kernel_position + 1));
  }
  // Second, I have to update merge cost value with previous kernel (if it exists).
  if(new_kernel_position == 0){
    return;
  }
  merge_cost_with_next_cluster_kernel_[new_kernel_position - 1] =
      CalculateDistanceBetweenClusterKernelAndTheirMerge(new_kernel_position - 1, new_kernel_position);
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
  UpdateMergeCostsListAfterMerging(index_of_ck_with_lowest_merge_cost_with_next_ck);
}

/** Updates merges costs list.
 * @brief Updates merges costs list.
 * @param merged_kernel_position
 */
void UnivariateListBasedClusterKernelAlgorithm::UpdateMergeCostsListAfterMerging(
    const unsigned int &merged_kernel_position) {
  // First I have to remove values for kernels that have been merged.
  auto iterator_of_merged_kernel = merge_cost_with_next_cluster_kernel_.begin() + merged_kernel_position;
  merge_cost_with_next_cluster_kernel_.erase(iterator_of_merged_kernel);
  merge_cost_with_next_cluster_kernel_.erase(iterator_of_merged_kernel);
  // Note that I could alternatively remove iterator_of_merged_kernel + 1 and iterator_of_merged_kernel,
  // but not in the reverse direction!
  // Now I have to do all the operations as if the only thing I did was adding merged kernel, thus...
  UpdateMergeCostsListAfterAddingKernel(merged_kernel_position);
}
