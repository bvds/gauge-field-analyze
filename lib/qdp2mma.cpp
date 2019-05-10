
/*
 * Created: 05-06-2017
 * Modified: Wed 07 Jun 2017 19:22:37 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 * ----------------------------------------------
 * Description:
 * Output Lattice fields in Mathematica format.
 * See 
 */

#include <qdp_converters.hpp>
#include <ostream>

namespace mma {
  
  template <typename Q>
  void write(ostream& to, const multi1d<Q>& a)
  {
    to << "{";
    for(int i=0; i<a.size; i++)
      {
	if(i>0)
	  to << ",";
	write(to, a[i]);
      }
    to << "}";
  }
  
  /* Color matrix element to Mathematica format */
  void write(ostream& to, QDP::ComplexD const &from)
  {
    to << QDP::real(from) << "+I*" << QDP::imag(from);
  }
  
  /* 
     See lib/qdp_converters.cpp from project [forked project].
     See QDP++ file lib/qdp_scalar_specific.cc, 
     especially the function writeArchiv 
  */



  /* Extracts the lattice gauge field and places it in a two dimensional array
   * where outer dim is direction and inner dim is coordinate
   */
  QDP::multi1d<QDP::multi1d<QDP::ColorMatrixD>>
  extract_colour_fields(QDP_Gauge_Field const &from)
  {
    QDP::multi1d<QDP::multi1d<QDP::ColorMatrixD>> to(Nd);
    
    for (int i = 0; i < Nd; ++i) {
      to[i].resize(QDP::Layout::sitesOnNode());
      QDP::QDP_extract(to[i], from[i], QDP::all);
    }
    
    return to;
  }


  /* 
     Lattice gauge field to Mathematica format
  */
  void write(ostream& to, const multi1d<LatticeColorMatrix>& u)
  {

    auto qdp_gauge_array = extract_colour_fields(from);

  // Might be worth doing a superindex instead, but it will require more
  // function calls and is in general less readable
  for (auto it = 0; it < L0; ++it)
    for (auto ix = 0; ix < L1; ++ix)
      for (auto iy = 0; iy < L2; ++iy)
        for (auto iz = 0; iz < L3; ++iz)
          for (auto mu = 0; mu < 4; ++mu)
            copy(qdp_gauge_array[mu][qdp_site_index(it, ix, iy, iz)],
                 to);

  void dump_gauge_lattice(, string& filename)
  {
  }

} // namespace mma
