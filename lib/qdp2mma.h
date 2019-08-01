
/*
 * Created: 05-06-2017
 * Modified: Wed 07 Jun 2017 16:07:44 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */

#ifndef QDP2MMA_H
#define QDP2MMA_H

#include <qdp.h>

namespace mma {

  void write(std::ostream& to, const QDP::ColorMatrix &from);
  void write(std::ostream& to, const QDP::Complex &from);
  void write(std::ostream& to, const QDP::Real &from);
  void write(std::ostream& to, const int &from);

  void dump_lattice(std::string& filename,
		    const QDP::multi1d<QDP::LatticeColorMatrix>& u,
		    const QDP::Real& beta);

} // namespace mma

#endif // QDP2MMA_H