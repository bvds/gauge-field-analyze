
/*
 * Created: 07-06-2017
 * Modified: Wed 07 Jun 2017 19:23:16 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */

#include <qdp_utilities.hpp>

namespace fastsum {

/* Computes the average plaquette of the QDP gauge config given. */
double average_plaquette(QDP_Gauge_Field const &gauge_field)
{
  QDP::Double plaq_sum{0.};

  for (auto mu = 0; mu < 4; ++mu)
    for (auto nu = 0; nu < mu; ++nu)
      plaq_sum += QDP::sum(QDP::real(QDP::trace(
          gauge_field[mu] * QDP::shift(gauge_field[nu], FORWARD, mu) *
          QDP::adj(QDP::shift(gauge_field[mu], FORWARD, nu)) *
          QDP::adj(gauge_field[nu]))));

  plaq_sum *= QDP::Double(1. / (QDP::Layout::vol() * 18.));

  return QDP::toDouble(plaq_sum);
}

} // namespace fastsum
