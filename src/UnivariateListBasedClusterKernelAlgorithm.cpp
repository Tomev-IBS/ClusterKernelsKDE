#include "UnivariateListBasedClusterKernelAlgorithm.h"
#include <cmath>

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
  UpdateMergeCostsList();
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
void UnivariateListBasedClusterKernelAlgorithm::UpdateMergeCostsList() {
  FillDomainForClusterKernelDistanceCalculation();
  merge_cost_with_next_cluster_kernel_.clear();
  for(auto i = 0; i < cluster_kernels_.size() - 1; ++i){
    merge_cost_with_next_cluster_kernel_.push_back(
      CalculateDistanceBetweenClusterKernelAndTheirMerge(i, i + 1)
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
  UpdateMergeCostsList();
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

/** Fills the domain for univariate case algorithm. We use the same approach as in case of DESDA algorithm.
 *  We take min and max values and then moves them by 5h (something like 5 sigma).
 * @brief Fills the domain for univariate case algorithm.
 */
void UnivariateListBasedClusterKernelAlgorithm::FillDomainForClusterKernelDistanceCalculation() {
  double minimal_value = FindMinimalValueOnDimension();
  double maximal_value = FindMaximalValueOnDimension();
  // Similarly to 5 sigma.
  minimal_value -= 5 * bandwidth_[0];
  maximal_value += 5 * bandwidth_[0];
  // We want hard 1000 points, so:
  double step_size = (maximal_value - minimal_value) / 1000;
  auto domain = std::vector<Point>();

  for(auto current_value = minimal_value; current_value < maximal_value; current_value += step_size){
    domain.push_back({current_value});
  }

  domain_for_cluster_kernel_distance_calculation_ = domain;
}

double UnivariateListBasedClusterKernelAlgorithm::FindMinimalValueOnDimension(const int &dimension) {
  if(cluster_kernels_.empty()){
    return 0;
  }

  double minimal_value_on_dimension = cluster_kernels_[0]->GetMean()[dimension];

  for(auto i = 1; i < cluster_kernels_.size(); ++i){
    if(cluster_kernels_[i]->GetMean()[dimension] < minimal_value_on_dimension){
      minimal_value_on_dimension = cluster_kernels_[i]->GetMean()[dimension];
    }
  }

  return minimal_value_on_dimension;
}

double UnivariateListBasedClusterKernelAlgorithm::FindMaximalValueOnDimension(const int &dimension) {
  if(cluster_kernels_.empty()){
    return 0;
  }

  double maximal_value_on_dimension = cluster_kernels_[0]->GetMean()[dimension];

  for(auto i = 1; i < cluster_kernels_.size(); ++i){
    if(cluster_kernels_[i]->GetMean()[dimension] > maximal_value_on_dimension){
      maximal_value_on_dimension = cluster_kernels_[i]->GetMean()[dimension];
    }
  }

  return maximal_value_on_dimension;
}

double UnivariateListBasedClusterKernelAlgorithm::CalculateDistanceBetweenClusterKernelAndTheirMerge(
    const int &first_ck_index, const int &second_ck_index) {
  auto merged_kernel = ClusterKernelPointer(cluster_kernels_[first_ck_index]->Merge(cluster_kernels_[second_ck_index].get()));
  merged_kernel->SetBandwidth(bandwidth_);

  double loss = 0;
  double domain_length = domain_for_cluster_kernel_distance_calculation_[domain_for_cluster_kernel_distance_calculation_.size() - 1][0] -
      domain_for_cluster_kernel_distance_calculation_[0][0];

  // Instead of integral, I'll use discrete sum over points.
  for(auto pt : domain_for_cluster_kernel_distance_calculation_){
    auto first_ck = cluster_kernels_[first_ck_index];
    auto second_ck = cluster_kernels_[second_ck_index];
    double current_loss = first_ck->GetCardinality() / bandwidth_[0] * first_ck->GetKernelValue({(pt[0] - first_ck->GetMean()[0]) / bandwidth_[0]})[0];
    current_loss += second_ck->GetCardinality() / bandwidth_[0] * second_ck->GetKernelValue({(pt[0] - second_ck->GetMean()[0]) / bandwidth_[0]})[0];
    current_loss -= merged_kernel->GetCardinality() / bandwidth_[0] * merged_kernel->GetKernelValue({(pt[0] - merged_kernel->GetMean()[0]) / bandwidth_[0]})[0];
    current_loss = pow(current_loss, 2);
    loss += current_loss;
  }
  loss /= domain_for_cluster_kernel_distance_calculation_.size();

  return sqrt(loss * domain_length);
}
