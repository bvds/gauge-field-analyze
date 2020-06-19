/*
 * Created: May 10, 2019
 * Author: Brett van de Sande (bvds@asu.edu)
 * ----------------------------------------------
 * Description:
 * Output Lattice fields in Mathematica format.
 */

#include "qdp2mma.h"
#include <iostream>
#include <fstream>

/*
 * The QDP++ macro CONFIG_LAYOUT is missing from the 
 * header files, so we have to brute-force it
 */
#if QDP_USE_LEXICO_LAYOUT == 1
#define CONFIG_LAYOUT "lexicographic"
#elif QDP_USE_CB2_LAYOUT == 1
#define CONFIG_LAYOUT "cb2"
#elif QDP_USE_CB3D_LAYOUT == 1
#define CONFIG_LAYOUT "cb3d"
#elif QDP_USE_CB32_LAYOUT == 1
#define CONFIG_LAYOUT "cb32"
#error "no appropriate layout defined"
#else

#endif

namespace mma {

  int count = 0;

  /* 
     See lib/qdp_converters.cpp from project qdp-to-openqcd on github.
     https://github.com/Irubataru/qdp-to-openqcd/blob/master/library/qdp_converters.cpp
     See QDP++ file lib/qdp_scalar_specific.cc, 
     especially the function writeArchiv 
  */
  template <typename Q>
  void write(std::ostream& to, const multi1d<Q>& a)
  {
    to << "{";
    for(int i=0; i<a.size(); i++)
      {
	if(i>0)
	  {
	    to << ",";
	    if(count>=4)
	      {
		to << std::endl;
		count = 0;
	      }
	  }
	write(to, a[i]);
      }
    to << "}";
  }

  // Extract color matrix components
  void write(std::ostream& to, const QDP::ColorMatrix &from)
  {
    // See QDP++/examples/reunit.cc
    QDP::multi1d<QDP::multi1d<QDP::Complex>> mm(Nc);

    for(int i=0; i<Nc; i++)
      {
	mm[i].resize(Nc);
	for(int j=0; j<Nc; j++)
	  mm[i][j] = QDP::peekColor(from, i, j);
      }
    write(to, mm);
  }

  /* Color matrix element to Mathematica format
     See QDP++/include/qdp_reality.h
  */
  void write(std::ostream& to, const QDP::Complex &from)
  {
    // Not sure why redundant type conversion is needed
    // Or why a cast to ordinary types are needed for camparision
    // operations.
    write(to, QDP::Real(QDP::real(from)));
    if(QDP::toWordType(QDP::Real(QDP::imag(from)))>=0.0)
      to  << "+";
    write(to, QDP::Real(QDP::imag(from)));
    to << "*I";
  }

  /* Float */
  void write(std::ostream& to, const QDP::Real &from)
  {
    count++;
    // Not sure why type conversion is needed here.
    if(QDP::toWordType(from)==0.0)
      to << from;
    else
      {
	int exp = floor(log10(fabs(QDP::toWordType(from))));
	// QDP operator != seems to be broken
	if(fabs(exp)>4)
	  to << from/pow(10.0, exp) << "*10^" << exp;
	else
	  to << from;
      }
  }

  /* Integer */
  void write(std::ostream& to, const int &from)
  {
    count++;
    to << from;
  }


  /* 
   *Extracts the lattice gauge field and places it in a two dimensional array
   * where outer dim is direction and inner dim is coordinate
   */
  QDP::multi1d<QDP::multi1d<QDP::ColorMatrix>>
  extract_colour_fields(const QDP::multi1d<QDP::LatticeColorMatrix> &from)
  {
    QDP::multi1d<QDP::multi1d<QDP::ColorMatrix>> to(Nd);

    for (int i = 0; i < Nd; ++i) {
      to[i].resize(QDP::Layout::sitesOnNode());
      /* QDP_extract(multi1d<OScalar<T> >& dest,
	 const OLattice<T>& src, const Subset& s)  */
      QDP::QDP_extract(to[i], from[i], QDP::all);
    }

    return to;
  }


  /* 
     Lattice gauge field to Mathematica format
  */
  void dump_lattice(std::string& filename,
		    const QDP::multi1d<QDP::LatticeColorMatrix>& u,
		    const Params& params)
  {
    std::ofstream to(filename);
    to.precision(15);  // Double, no matter what
    to << std::fixed;  // Fixed point:  handle exponents manually

    // A lot of this is also in the XML file.
    to << "beta=" << params.beta << "; ";
    to << "nd=" << Nd << "; nc=" << Nc << ";" << std::endl;
    to << "latticeDimensions=";
    write(to, QDP::Layout::lattSize());
    to << ";" << std::endl;
    to << "latticeLayout=\"" << CONFIG_LAYOUT << "\";" << std::endl;

    to << "latticeBC=\"" << params.GaugeBC << "\"; " << std::endl;
    if(params.GaugeBC == "TRANS_GAUGEBC") {
        count = 0;
        to << "transverseBCPhases=";
        write(to, params.phases);
        to << "; " << std::endl;
        count = 0;
        to << "transverseBCDirections=";
        // Shift to match Mathematica array numbering.
        write(to, params.transverseDirs+1);
        to << ";" << std::endl;
    } else {
        to << "transverseBCPhases=Null; ";
        to << "transverseBCDirections=Null;" << std::endl;
    }

    to << "gaugeField=";

    auto qdp_gauge_array = extract_colour_fields(u);
    write(to, qdp_gauge_array);

    to << ";" << std::endl;

    /*
     * Lattice Coordinates
     * Instead of figuring out how checkerboard works,
     * just provide a lookup table.
     */
    to << "Clear[linearSiteIndex];" << std::endl;
    for(int linear=0; linear<QDP::Layout::vol(); ++linear)
      {
	// Get the true lattice coord of this linear site index
	multi1d<int> coord = Layout::siteCoords(0, linear);
	// Shift coordinates by 1 to match Mathematica array numbering
	for(int i=0; i<Nd; i++)
	  coord[i]++;
	count=0;
	to << "linearSiteIndex[";
	write(to, coord);
	to << "]=" << linear+1 << ";" << std::endl;
      }

    to.close();
  }

} // namespace mma
