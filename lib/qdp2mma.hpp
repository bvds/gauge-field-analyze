
/*
 * Created: 05-06-2017
 * Modified: Wed 07 Jun 2017 16:07:44 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */

#ifndef QDP_CONVERTERS_HPP
#define QDP_CONVERTERS_HPP

#include <qdp.h>
#include <types.hpp>

namespace fastsum {

void copy(QDP::ColorMatrixD const& from);
void copy(QDP_Gauge_Field const& from);

} // namespace fastsum 

#endif /* QDP_CONVERTERS_HPP */
