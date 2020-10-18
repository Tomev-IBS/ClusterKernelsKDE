#include "WeightedUnivariateListBasedClusterKernelAlgorithm.h"
#include <iostream>

namespace WeightedUnivariateListBasedClusterKernelAlgorithmUtilities{
  /** Multiplies the given vector by the given scalar. Vector here is defined by a point.
  * @brief Multiplies the given vector by the given scalar.
  * @param point - A point that defines a vector that will be multiplied.
  * @param scalar - A scalar by which vector will be multiplied.
  * @return Multiplied vector.
  */
  Point MultiplyVectorByScalar(const Point &point, const double &scalar){
    Point multiplied_vector = {};

    for(auto value : point){
      multiplied_vector.push_back(scalar * value);
    }

    return multiplied_vector;
  }

  /** Return sum of two vectors given by points (vectors of doubles). For more physical view
 *  one can assume that the vectors have one end in given point and the second at 0.
 *  @brief Return sum of two vectors given by points (vectors of doubles).
 *  @param first_point - A point defining the first vector.
 *  @param second_point - A point defining the second vector.
 *  @return Sum of the vectors.
 */
  Point SumVectors(const Point &first_point, const Point &second_point){
    Point sum = {};
    // Check if dimensions are right. Return empty, if not.
    if(second_point.size() != first_point.size()){
      return sum;
    }

    for(auto i = 0; i < first_point.size(); ++i){
      sum.push_back(first_point[i] + second_point[i]);
    }

    return sum;
  }
}

using namespace WeightedUnivariateListBasedClusterKernelAlgorithmUtilities;

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
  if(!cluster_kernels_.empty()) {
    new_cluster_kernel->RescaleWeight(weight_modifier_);
  }
  // This is for univariate case. The idea is to place new cluster kernel in between
  // means closest to it.
  double new_cluster_kernel_value = new_cluster_kernel->GetMean()[0];
  unsigned int new_kernel_position = 0;
  for(auto kernel : cluster_kernels_){
    if(kernel->GetMean()[0] > new_cluster_kernel_value) {
      break;
    }
    ++new_kernel_position;
  }
  new_cluster_kernel->SetBandwidth(bandwidth_);
  cluster_kernels_.insert(cluster_kernels_.begin() + new_kernel_position, new_cluster_kernel);

  /*
  std::cout << "Current kernels means:\n";
  for(auto kernel : cluster_kernels_){
    std::cout << "\t" << kernel->GetMean()[0] << std::endl;
  }
  */
}

/** Calculates value of KDE in given point. Note, that due to using cluster kernels this is
 *  an approximation (in comparison to standard KDE).
 * @brief Calculates value of KDE in given point.
 * @param x - point for which value is calculated.
 * @return Value of KDE in given point.
 */
Point WeightedUnivariateListBasedClusterKernelAlgorithm::GetValue(const Point &x) {
  if(cluster_kernels_.empty() || number_of_parsed_elements_ == 0){
    return {};
  }

  Point value_at_point = MultiplyVectorByScalar(cluster_kernels_[0]->GetValue(x), cluster_kernels_[0]->GetWeight());

  for(auto i = 1; i < cluster_kernels_.size(); ++i){
    auto addend_from_kernel = MultiplyVectorByScalar(cluster_kernels_[i]->GetValue(x), cluster_kernels_[i]->GetWeight());
    value_at_point = SumVectors(value_at_point, addend_from_kernel);
  }

  return value_at_point;
}
