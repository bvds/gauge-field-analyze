
/*
 * Created: 05-06-2017
 * Modified: Wed 07 Jun 2017 19:22:37 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 * ----------------------------------------------
 * Description:
 * Functions that converts between the QDP and the OpenQCD storage format.
 * Heavily inspired by the Mainz "observer" code.
 */

#include <qdp_converters.hpp>

namespace fastsum {

namespace {

/* Copy the complex variable of a single colour index.
 * Used instead of a loop because OpenQCD does not store the colour matrix as an
 * array.
 */
void copy_colour_index(QDP::ColorMatrixD const &from, complex_dble &to, int i,
                       int j)
{
  static QDP::ComplexD tmp;

  tmp = QDP::peekColor(from, i, j);
  to.re = QDP::toDouble(QDP::real(tmp));
  to.im = QDP::toDouble(QDP::imag(tmp));
}

/* Extracts the lattice gauge field and places it in a two dimensional array
 * where outer dim is direction and inner dim is 4D coordinate
 */
QDP::multi1d<QDP::multi1d<QDP::ColorMatrixD>>
extract_colour_fields(QDP_Gauge_Field const &from)
{
  QDP::multi1d<QDP::multi1d<QDP::ColorMatrixD>> to(4);

  for (auto i = 0ul; i < 4; ++i) {
    to[i].resize(QDP::Layout::sitesOnNode());
    QDP::QDP_extract(to[i], from[i], QDP::all);
  }

  return to;
}

/* The superindex of pos (t,x,y,z) for the QDP storage */
int qdp_site_index(int t, int x, int y, int z)
{
  auto index_array = QDP::multi1d<int>(4);

  index_array[0] = x;
  index_array[1] = y;
  index_array[2] = z;
  index_array[3] = t;

  return QDP::Layout::linearSiteIndex(index_array);
}

/* The superindex of pos (t,x,y,z) for the OpenQCD storage */
int openqcd_site_index(int t, int x, int y, int z)
{
  return ipt[z + L3 * (y + L2 * (x + L1 * t))];
}

/* Checks if a site is odd or even */
bool is_odd_site(int t, int x, int y, int z)
{
  return ((t + x + y + z) % 2 == 1);
}

/* Returns the index of the gauge link at pos (t,x,y,z) in dir mu in the openqcd
 * storage format. It assumes that the program runs in serial and that all
 * openqcd geometry structures exist. The mu index is assumed to be of the
 * (x,y,z,t) order of QDP and must therefore be converted.
 */
int openqcd_gauge_index(int t, int x, int y, int z, int mu)
{
  auto site_index = openqcd_site_index(t, x, y, z);
  int omu = mu + 1 - 4*(mu == 3);

  if (is_odd_site(t, x, y, z))
    return 8 * (site_index - VOLUME / 2) + 2 * omu;
  else
    return 8 * (iup[site_index][omu] - VOLUME / 2) + 2 * omu + 1;
}

}

/* Single colour matrix copy function */
void copy(QDP::ColorMatrixD const &from, su3_dble &to)
{
  QDP::ComplexD tmp;

  copy_colour_index(from, to.c11, 0, 0);
  copy_colour_index(from, to.c12, 0, 1);
  copy_colour_index(from, to.c13, 0, 2);

  copy_colour_index(from, to.c21, 1, 0);
  copy_colour_index(from, to.c22, 1, 1);
  copy_colour_index(from, to.c23, 1, 2);

  copy_colour_index(from, to.c31, 2, 0);
  copy_colour_index(from, to.c32, 2, 1);
  copy_colour_index(from, to.c33, 2, 2);
}

/* Lattice gauge field copy function */
void copy(QDP_Gauge_Field const &from, OpenQCD_Gauge_Field &to)
{
  if (!OpenQCD::Is_Initialised())
    throw std::runtime_error{
        "Trying to copy fields without initialising OpenQCD"};

  auto qdp_gauge_array = extract_colour_fields(from);

  // Might be worth doing a superindex instead, but it will require more
  // function calls and is in general less readable
  for (auto it = 0; it < L0; ++it)
    for (auto ix = 0; ix < L1; ++ix)
      for (auto iy = 0; iy < L2; ++iy)
        for (auto iz = 0; iz < L3; ++iz)
          for (auto mu = 0; mu < 4; ++mu)
            copy(qdp_gauge_array[mu][qdp_site_index(it, ix, iy, iz)],
                 to[openqcd_gauge_index(it, ix, iy, iz, mu)]);
}

} // namespace fastsum
