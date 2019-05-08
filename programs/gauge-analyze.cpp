
/*
 * Created: 02-06-2017
 * Modified: Wed 07 Jun 2017 16:34:22 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */

// Needed to define all OpenQCD global arrays
#define MAIN_PROGRAM

#include "gauge-analyze.hpp"
#include <iomanip>

using namespace fastsum;

const std::string xml_geometry_tag = "/HMC/ProgramInfo/Setgeom";
const std::string xml_lattice_size_tag = "latt_size";

int main(int argc, char **argv)
{
  QDP_Instance_Wrapper qdp_instance{argc, argv};

  if (QDP_Instance_Wrapper::MPI::comm_size() > 1 or NPROC != 1) {
    std::cerr << "This program only works in serial" << std::endl;
    return 1;
  }

  Program_Parameters params;

  try {
    params = parse_input_arguments(argc, argv);
  } catch (std::exception &err) {
    std::cerr << "Error: " << err.what() << "\n" << std::endl;
    print_help_message();
    return 1;
  }

  OpenQCD::Initialise();
  init_qdp_lattice_geometry();

  QDP_Gauge_Field qdp_gauge_field(4);

  try {
    read_qdp_gauge_field(qdp_gauge_field, params.input_file);
  } catch (std::exception &err) {
    std::cerr << err.what() << std::endl;
    return 1;
  }

  std::cout << "QDP config imported, average plaquette:     "
            << std::setprecision(16) << average_plaquette(qdp_gauge_field)
            << std::endl;

  OpenQCD_Gauge_Field openqcd_gauge_field{udfld()};

  try {
    copy(qdp_gauge_field, openqcd_gauge_field);
    export_openqcd_config(params.output_file);
  } catch (std::exception &err) {
    std::cerr << err.what() << std::endl;
    return 1;
  }

  std::cout << "OpenQCD config exported, average plaquette: "
            << std::setprecision(16) << average_plaquette(openqcd_gauge_field)
            << std::endl;
}

namespace fastsum {

void init_qdp_lattice_geometry()
{
  QDP::multi1d<int> nrow(4);
  nrow[0] = NPROC1 * L1;
  nrow[1] = NPROC2 * L2;
  nrow[2] = NPROC3 * L3;
  nrow[3] = NPROC0 * L0;

  QDP::Layout::setLattSize(nrow);
  QDP::Layout::create();
}

void read_qdp_gauge_field(QDP_Gauge_Field &gauge_field, std::string filename)
{
  QDP::XMLReader lime_xml_header{};
  QDP::QDPFileReader lime_reader{lime_xml_header, filename, QDP::QDPIO_SERIAL};

  check_input_geometry(lime_xml_header);

  QDP::XMLReader lime_xml_record{};
  QDP::read(lime_reader, lime_xml_record, gauge_field);

  if (lime_reader.bad())
    throw std::runtime_error{"Errors reading input config file"};
}

void check_input_geometry(QDP::XMLReader &lime_xml_header)
{
  QDP::XMLReader lime_geom{lime_xml_header, xml_geometry_tag};

  QDP::multi1d<int> lime_rows;
  QDP::read(lime_geom, xml_lattice_size_tag, lime_rows);

  if ((lime_rows[0] != NPROC1 * L1) or (lime_rows[1] != NPROC2 * L2) or
      (lime_rows[2] != NPROC3 * L3) or (lime_rows[3] != NPROC0 * L0))
    throw std::runtime_error{"openqcd and input file geometry mismatch"};
}

void print_help_message()
{
  std::cout << "Usage:\n"
            << "qdp_to_openqcd [input_filename] [output_filename]" << std::endl;
}

Program_Parameters parse_input_arguments(int argc, char **argv)
{
  if (argc < 3)
    throw std::runtime_error{"Not enough arguments"};

  auto result = Program_Parameters{{argv[1]}, {argv[2]}};

  std::ifstream ifs{result.input_file};
  if (!ifs)
    throw std::runtime_error{"Cannot open \"" + result.input_file +
                             "\" for reading"};
  ifs.close();

  std::ofstream ofs{result.output_file};
  if (!ofs)
    throw std::runtime_error{"Cannot open \"" + result.output_file +
                             "\" for writing"};
  ofs.close();

  return result;
}

} // namespace fastsum
