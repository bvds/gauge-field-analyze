
/*
 * Created: 07-06-2017
 * Modified: Wed 07 Jun 2017 19:21:36 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */

#include <openqcd_utilities.hpp>

extern "C" {
#include "uflds.h"
}

namespace fastsum {

/* Compute the average plaquette of the OpenQCD gauge field. It takes it as an
 * argument, but it uses the global field as stored in udfld(). The argument is
 * simply there for function overloading reasons.
 */
double average_plaquette(OpenQCD_Gauge_Field const &)
{
  static constexpr double normalisation = 18. * VOLUME;
  if (!OpenQCD::Is_Initialised())
    throw std::runtime_error{
        "Cannot measure the plaquette before OpenQCD is initialised"};

  return plaq_wsum_dble(0) / normalisation;
}

} // namespace fastsum
