
/*
 * Created: 02-06-2017
 * Modified: Wed 07 Jun 2017 16:22:59 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */

#ifndef GAUGE_ANALYZE_H
#define GAUGE_ANALYZE_H

#include <iostream>
#include <qdp.h>
#include "gauge_io.h"
#include "qdp2mma.h"
#include "qdp_utilities.h"

namespace fastsum {

struct Program_Parameters
{
  std::string input_file;
  std::string output_file;
};

void init_qdp_lattice_geometry();
void read_qdp_gauge_field(QDP::multi1d<QDP::LatticeColorMatrix> &gauge_field,
			  std::string filename);
void check_input_geometry(QDP::XMLReader &lime_xml_header);
void print_help_message();
Program_Parameters parse_input_arguments(int arc, char **argv);

} // namespace fastsum

#endif /* GAUGE_ANALYZE_H */
