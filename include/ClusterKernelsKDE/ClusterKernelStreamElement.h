#ifndef CLUSTERKERNELSKDE_CLUSTERKERNELSTREAMELEMENT_H
#define CLUSTERKERNELSKDE_CLUSTERKERNELSTREAMELEMENT_H

/** This is an abstract class for stream elements for Cluster Kernels stream elements. It contains every functionality
 * that an stream element has to provide in order for Cluster Kernels algorithm to be able to properly use it.
 * @brief This is an abstract class for stream elements for Cluster Kernels stream elements.
 */

#include <vector>

typedef std::vector<double> point;

class ClusterKernelStreamElement{
  public:
    virtual point GetMean() = 0;
};

#endif //CLUSTERKERNELSKDE_CLUSTERKERNELSTREAMELEMENT_H
