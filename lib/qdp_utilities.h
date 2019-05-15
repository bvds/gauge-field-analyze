
/*
 * Created: 07-06-2017
 * Modified: Wed 07 Jun 2017 16:08:23 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */

#ifndef QDP_UTILITIES_H
#define QDP_UTILITIES_H

#include <qdp.h>

namespace fastsum {

double average_plaquette(const QDP::multi1d<QDP::LatticeColorMatrix>& u);

} // namespace fastsum 


#endif /* QDP_UTILITIES_H */
