#include "ClusterKernelsAlgorithm.h"

#include <cmath>
#include <algorithm>

bool AlmostEqual(const double &a, const double &b) {
  // I'm using method from http://realtimecollisiondetection.net/blog/?p=89.
  double absolute_tolerance = 1e-5;
  double relative_tolerance = 1e-5;

  return fabs(a - b) <= std::max(absolute_tolerance,
                                    relative_tolerance * std::max(fabs(a), fabs(b)));
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

/** Return difference of two vectors given by points (vectors of doubles). For more physical view
 *  one can assume that the vectors have one end in given point and the second at 0.
 *  @brief Return difference of two vectors given by points (vectors of doubles).
 *  @param first_point - A point defining the first vector.
 *  @param second_point - A point defining the second vector.
 *  @return Difference of the vectors.
 */
Point SubtractVectors(const Point &first_point, const Point &second_point){
  Point difference = {};
  // Check if dimensions are right. Return empty, if not.
  if(second_point.size() != first_point.size()){
    return difference;
  }

  for(auto i = 0; i < first_point.size(); ++i){
    difference.push_back(first_point[i] + second_point[i]);
  }

  return difference;
}

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

/** Calculates length of the vector. This method in particular calculates euclidean distance.
 * @brief Calculates length of the vector.
 * @param pt - Point defining the vector.
 * @return - Given vectors length.
 */
double GetVectorsLength(const Point &pt){
  double vectors_length_squared = 0;
  for(auto value : pt){
    vectors_length_squared += pow(value, 2);
  }
  return sqrt(vectors_length_squared);
}

ClusterKernelsAlgorithm::ClusterKernelsAlgorithm(const int &m,
                                                 ClusterKernel*(*cluster_kernel_factory_method)(ClusterKernelStreamElement *))
  : maximal_number_of_cluster_kernels_(m),
  cluster_kernel_factory_method_(cluster_kernel_factory_method)
{}

void ClusterKernelsAlgorithm::PerformStep(ClusterKernelStreamElement *stream_element) {
  UpdateBandwidth(stream_element);
  /*
   * The algorithm usually starts with finding out if there exist a cluster kernel with
   * mean equal to stream elements mean. Here I begin with finding which index is it or
   * in case there are bit "same-mean" cluster kernels return the size of cluster kernels
   * index (would-be next index).
    */

  int index_of_cluster_kernel_with_same_mean_as_stream_element =
      FindIndexOfClusterKernelWithSameMeanAsStreamElement(stream_element);
  // Decide if new cluster kernel should be added or cluster kernel should be updated.
  if(index_of_cluster_kernel_with_same_mean_as_stream_element == cluster_kernels_.size()){
    AddNewClusterKernel(stream_element);
  }
  else{
    cluster_kernels_[index_of_cluster_kernel_with_same_mean_as_stream_element]->Update(stream_element);
  }

  // Merge cluster kernels if there are too many
  if(cluster_kernels_.size() > maximal_number_of_cluster_kernels_) {
    MergeClusterKernelsWithTheLowestMergeCost();
  }
}

/** This method finds and returns a Cluster Kernel with the "same" mean as the stream element. It
 * uses implemented euclidean distance. Note, that if no Cluster Kernel had the "same" mean
 * as the stream element, the function will return the size of cluster kernels vector.
 *
 * @brief This method finds and returns a Cluster Kernel closest to given stream element.
 * @param streamElement
 * @return
 */
int ClusterKernelsAlgorithm::FindIndexOfClusterKernelWithSameMeanAsStreamElement(ClusterKernelStreamElement *streamElement) const {
  int index = 0;

  for(;index < cluster_kernels_.size(); ++index){
    if(AlmostEqual(0, CalculateDistanceBetweenPoints(streamElement->GetMean(),
                                                     cluster_kernels_[index]->GetMean()))){
      return index;
    }
  }

  return index;
}

/** Method for calculating distance between points. In this case euclidean distance is calculated
 * as it's one of the easiest to understand and the most well known.
 *
 * TR TODO: Consider having distance calculating object in case different distance measures
 * would be required.
 *
 * @brief Method for calculating distance between points.
 * @param point1
 * @param point2
 * @return Distance between two given points.
 */
double ClusterKernelsAlgorithm::CalculateDistanceBetweenPoints(Point point1, Point point2) const {
  // Check if dimensions match.
  if(point1.size() != point2.size()){
    return -1; // Return -1 if not (which clearly cannot be distance thus indicates something's up)
  }
  double distance = 0;

  for(auto i = 0; i < point1.size(); ++i){
    distance += pow(point1[i] - point2[i], 2);
  }

  return sqrt(distance);
}

/** Adds new cluster kernel to cluster kernels vector. This is done using cluster kernel factory
 * method provided in the constructor of this object.
 * @brief Adds new cluster kernel to cluster kernels vector.
 * @param stream_element Stream element basing on which new cluster kernel is constructed.
 */
void ClusterKernelsAlgorithm::AddNewClusterKernel(ClusterKernelStreamElement *stream_element) {
  cluster_kernels_.push_back(
    ClusterKernelPointer(cluster_kernel_factory_method_(stream_element))
  );
}

/** Merges two cluster kernels with the lowest merge cost. Instead of calculating l2 cost, as
 *  proposed in original work, I'll use total variance distance in order to determine which merged cluster kernel
 *  is the closest to it's part. All the xs used for this calculation will basically be means of cluster
 *  kernels.
 *
 * @brief Merges two cluster kernels with the lowest merge cost.
 */
void ClusterKernelsAlgorithm::MergeClusterKernelsWithTheLowestMergeCost() {
  FillDomainForClusterKernelDistanceCalculation();
  // Finding i and j such that i-th and j-th cluster kernel has the lowest distance from their merge.
  double minimal_distance = 1; // Initialize with maximum TV Distance value.
  unsigned int first_cluster_kernel_index = 0;
  unsigned int second_cluster_kernel_index = 0;
  for(auto i = 0; i < cluster_kernels_.size(); ++i){
    for(auto j = i + 1; j < cluster_kernels_.size(); ++j){
      double current_distance = CalculateDistanceBetweenClusterKernelAndTheirMerge(i, j);
      if(current_distance < minimal_distance){
        minimal_distance = current_distance;
        first_cluster_kernel_index = i;
        second_cluster_kernel_index = j;
      }
    }
  }
  MergeClusterKernels(first_cluster_kernel_index, second_cluster_kernel_index);
}

/** This method fills the domain for cluster kernel distance calculation with points. Although it may be
 * too little, in this case I'll use only means of cluster kernels, as a valid (and most important!)
 * points of the domain.
 *
 * @brief This method fills the domain for cluster kernel distance calculation with points.
 */
void ClusterKernelsAlgorithm::FillDomainForClusterKernelDistanceCalculation() {
  domain_for_cluster_kernel_distance_calculation_.clear();
  for(auto cluster_kernel : cluster_kernels_){
    domain_for_cluster_kernel_distance_calculation_.push_back(cluster_kernel->GetMean());
  }
}

/** Calculates distance between merged cluster kernel and it's parts. Note, that instead of
 *  loss function presented in the original work I use total variance distance.
 * @brief Calculates distance between merged cluster kernel and it's parts.
 * @param first_ck_index - Index of first cluster kernel.
 * @param second_ck_index - Index of second cluster kernel.
 * @return TV Distance between cluster kernels and their merge.
 */
double ClusterKernelsAlgorithm::CalculateDistanceBetweenClusterKernelAndTheirMerge(const int &first_ck_index,
                                                                                   const int &second_ck_index) {
  auto merged_kernel = std::shared_ptr<ClusterKernel>(cluster_kernels_[first_ck_index]->Merge(*cluster_kernels_[second_ck_index]));
  // Initialize the point denoting difference.
  Point current_difference = SumVectors(
      cluster_kernels_[first_ck_index]->GetValue(domain_for_cluster_kernel_distance_calculation_[0]),
      cluster_kernels_[second_ck_index]->GetValue(domain_for_cluster_kernel_distance_calculation_[0])
      );

  current_difference = SubtractVectors(current_difference, merged_kernel->GetValue(domain_for_cluster_kernel_distance_calculation_[0]));

  double difference = GetVectorsLength(current_difference);

  for(auto i = 1; i < cluster_kernels_.size(); ++i){
    current_difference = SumVectors(
        cluster_kernels_[first_ck_index]->GetValue(domain_for_cluster_kernel_distance_calculation_[i]),
        cluster_kernels_[second_ck_index]->GetValue(domain_for_cluster_kernel_distance_calculation_[i])
                                   );
    current_difference = SubtractVectors(current_difference, merged_kernel->GetValue(domain_for_cluster_kernel_distance_calculation_[i]));
    difference += GetVectorsLength(current_difference);
  }

  return 0.5 * difference;
}

/** Merges i-th and j-th cluster kernel. Algorithm is set in such way, that i < j, thus
 *  merged kernel should replace i-th kernel and j-th kernel should be removed afterwards.
 *  @brief Merges i-th and j-th cluster kernel.
 *  @param first_kernel_index - Index of the first cluster kernel to merge.
 *  @param second_kernel_index - Index of the second cluster kernel to merge.
 */
void ClusterKernelsAlgorithm::MergeClusterKernels(const unsigned int &first_kernel_index,
                                                  const unsigned int &second_kernel_index) {
  cluster_kernels_[first_kernel_index] = std::shared_ptr<ClusterKernel>(
      cluster_kernels_[first_kernel_index]->Merge(*cluster_kernels_[second_kernel_index].get()));
  cluster_kernels_.erase(cluster_kernels_.begin() + second_kernel_index);
}

/**
 * @brief This methods updates bandwidth according to 2008 Cluster Kernels work.
 * @param stream_element - Actually parsed stream element.
 */
void ClusterKernelsAlgorithm::UpdateBandwidth(ClusterKernelStreamElement *stream_element) {
  UpdateVariationEstimator(stream_element);

  bandwidth_.clear();
  if(number_of_parsed_elements_ == 1){
    // In the first step bandwidth should be initialized with 1.
    for(auto value : variation_estimator_) {
      bandwidth_.push_back(1);
    }
  }
  else {
    for(auto value : variation_estimator_) {
      bandwidth_.push_back(
          bandwidth_coefficient_ * sqrt(value) * pow(number_of_parsed_elements_, elements_number_coefficient));
    }
  }
}

/** Update variation using Welford's online formula. As a side effect and a required step it
 *  also updates mean estimator.
 *  @brief Update variation using Welford's online formula.
 *  @param stream_element - Actually parsed stream element.
 */
void ClusterKernelsAlgorithm::UpdateVariationEstimator(ClusterKernelStreamElement *stream_element) {
  ++number_of_parsed_elements_;
  if(number_of_parsed_elements_ == 0){
    last_step_elements_mean_ = stream_element->GetMean();
    for(auto dimension = 0; dimension < last_step_elements_mean_.size(); ++dimension){
      variation_estimator_.push_back(0); // Initialize variation estimator.
    }
  }
  else{
    Point elements_mean = stream_element->GetMean();
    Point current_mean = {};
    Point current_variation = {};
    for(auto dimension = 0; dimension < last_step_elements_mean_.size(); ++dimension){
      // This loop is done for each dimension separately. Dimension values are therefore
      // meant with respect to currently parsed dimension.
      double dimensions_mean = last_step_elements_mean_[dimension];
      dimensions_mean += (elements_mean[dimension] - last_step_elements_mean_[dimension]) / number_of_parsed_elements_;
      current_mean.push_back(dimensions_mean);
      double dimensions_variation = variation_estimator_[dimension];
      dimensions_variation += ((elements_mean[dimension] - last_step_elements_mean_[dimension])
          * (elements_mean[dimension] - dimensions_mean) - last_step_elements_mean_[dimension]) / number_of_parsed_elements_;
      current_variation.push_back(dimensions_variation);
    }
    last_step_elements_mean_ = current_mean;
    variation_estimator_ = current_variation;
  }
}

/** Calculates value of KDE in given point. Note, that due to using cluster kernels this is
 *  an approximation (in comparison to standard KDE).
 * @brief Calculates value of KDE in given point.
 * @param x - point for which value is calculated.
 * @return Value of KDE in given point.
 */
Point ClusterKernelsAlgorithm::GetValue(const Point &x) {
  if(cluster_kernels_.empty() || number_of_parsed_elements_ == 0){
    return {};
  }

  Point value_at_point = cluster_kernels_[0]->GetValue(x);

  for(auto i = 1; i < cluster_kernels_.size(); ++i){
    value_at_point = SumVectors(value_at_point, cluster_kernels_[i]->GetValue(x));
  }

  for(auto i = 0; i < value_at_point.size(); ++i){
    value_at_point[i] /= number_of_parsed_elements_;
  }

  return value_at_point;
}

Point ClusterKernelsAlgorithm::CalculateMergedClusterKernelMean(const unsigned int &first_kernel_index,
                                                                const unsigned int &second_kernel_index) {
  double weights_sum = cluster_kernels_[first_kernel_index]->GetWeight() +
      cluster_kernels_[second_kernel_index]->GetWeight();

  // In case there's problem with weights.
  if(AlmostEqual(weights_sum, 0)){
    return {};
  }

  Point first_weighted_mean = MultiplyVectorByScalar(
      cluster_kernels_[first_kernel_index]->GetMean(),
      cluster_kernels_[first_kernel_index]->GetWeight()
                                                    );
  Point second_weighted_mean = MultiplyVectorByScalar(
      cluster_kernels_[second_kernel_index]->GetMean(),
      cluster_kernels_[second_kernel_index]->GetWeight()
                                                     );

  return MultiplyVectorByScalar(SumVectors(first_weighted_mean, second_weighted_mean), 1.0 / weights_sum);
}
