
/*
 * Created: 05-06-2017
 * Modified: Wed 07 Jun 2017 16:07:44 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */

#ifndef QDP_CONVERTERS_HPP
#define QDP_CONVERTERS_HPP

#include <qdp.h>
#include <types.hpp>
#include <openqcd_utilities.hpp>

namespace fastsum {

void copy(QDP::ColorMatrixD const&from, su3_dble &to);
void copy(QDP_Gauge_Field const& from, OpenQCD_Gauge_Field &to);

} // namespace fastsum 

#endif /* QDP_CONVERTERS_HPP */
