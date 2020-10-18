#ifndef CLUSTERKERNELSKDE_CLUSTERKERNELSALGORITHM_H
#define CLUSTERKERNELSKDE_CLUSTERKERNELSALGORITHM_H

#include <vector>

#include "ClusterKernelStreamElement.h"
#include "ClusterKernel.h"


class ClusterKernelsAlgorithm : public RealValuedFunction{

  public:
    ClusterKernelsAlgorithm(const int &m,
                            ClusterKernel*(*cluster_kernel_factory_method)(ClusterKernelStreamElement *stream_element));
    virtual void PerformStep(ClusterKernelStreamElement *stream_element);
    Point GetValue(const Point &pt) override;
  protected:
    int maximal_number_of_cluster_kernels_ = 0;
    std::vector<ClusterKernelPointer> cluster_kernels_;
    ClusterKernel* (*cluster_kernel_factory_method_)(ClusterKernelStreamElement *stream_element);
    std::vector<Point> domain_for_cluster_kernel_distance_calculation_;

    // Bandwidth related variables
    std::vector<double> bandwidth_ = {};
    unsigned long long number_of_parsed_elements_ = 0;
    Point last_step_elements_mean_ = {};
    Point variation_estimator_ = {};
    // "Magic values" from the paper
    double bandwidth_coefficient_ = 1.06;
    double elements_number_coefficient = -0.2;

    void UpdateBandwidth(ClusterKernelStreamElement *stream_element);
    void UpdateVariationEstimator(ClusterKernelStreamElement *stream_element);
    int FindIndexOfClusterKernelWithSameMeanAsStreamElement(ClusterKernelStreamElement *stream_element) const;
    double CalculateDistanceBetweenPoints(Point point1, Point point2) const;
    virtual void AddNewClusterKernel(ClusterKernelStreamElement *stream_element);
    virtual void MergeClusterKernelsWithTheLowestMergeCost();
    virtual void FillDomainForClusterKernelDistanceCalculation();
    virtual double CalculateDistanceBetweenClusterKernelAndTheirMerge(const int &first_ck_index, const int &second_ck_index);
    void MergeClusterKernels(const unsigned int &first_kernel_index, const unsigned int &second_kernel_index);
};

#endif //CLUSTERKERNELSKDE_CLUSTERKERNELSALGORITHM_H
