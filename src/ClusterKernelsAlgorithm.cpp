#include "ClusterKernelsAlgorithm.h"

#include <cmath>
#include <algorithm>

bool AlmostEqual(const double &a, const double &b) {
  // I'm using method from http://realtimecollisiondetection.net/blog/?p=89.
  double absolute_tolerance = 1e-5;
  double relative_tolerance = 1e-5;

  return abs(a - b) <= std::max(absolute_tolerance,
                                    relative_tolerance * std::max(abs(a), abs(b)));
}

ClusterKernelsAlgorithm::ClusterKernelsAlgorithm(const int &m, ClusterKernel*(*cluster_kernel_factory_method)(ClusterKernelStreamElement *stream_element))
  : maximal_number_of_cluster_kernels_(m),
  cluster_kernel_factory_method_(cluster_kernel_factory_method)
{}

void ClusterKernelsAlgorithm::PerformStep(ClusterKernelStreamElement *stream_element) {
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
double ClusterKernelsAlgorithm::CalculateDistanceBetweenPoints(point point1, point point2) const {
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
 *  is the closest to it's part.
 *
 * @brief Merges two cluster kernels with the lowest merge cost.
 */
void ClusterKernelsAlgorithm::MergeClusterKernelsWithTheLowestMergeCost() {
  CalculateDomainForClusterKernelCalculation();
}

void ClusterKernelsAlgorithm::CalculateDomainForClusterKernelCalculation() {
  domain_for_cluster_kernel_distance_calculation_.clear();
  // First find min and max mean in within cluster kernels.

}
