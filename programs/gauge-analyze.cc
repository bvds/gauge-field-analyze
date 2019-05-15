
/*
 * Created: 02-06-2017
 * Modified: Wed 07 Jun 2017 16:34:22 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */


#include "gauge-analyze.h"
#include <iomanip>

struct Program_Parameters
{
  std::string input_file;
  std::string output_file;
};

Program_Parameters parse_input_arguments(int argc, char **argv)
{
  if(argc < 3)
    throw std::runtime_error{"Not enough arguments"};

  auto result = Program_Parameters {{argv[1]}, {argv[2]}};

  std::ifstream ifs {result.input_file};
  if (!ifs)
    throw std::runtime_error {"Cannot open \"" + result.input_file + "\" for reading"};
  ifs.close();

  std::ofstream ofs {result.output_file};
  if (!ofs)
    throw std::runtime_error {"Cannot open \"" + result.output_file + "\" for writing"};
  ofs.close();

  return result;
}


int main(int argc, char **argv)
{
  // This only works for scalar!

  // Taken from Chroma file lib/init/chroma_init.cc
#if defined QDPJIT_IS_QDPJITPTX || defined QDPJIT_IS_QDPJITNVVM
  if (! QDP::QDP_isInitialized())
    QDP::QDP_initialize_CUDA(&argc, &argv);
#else
  if (! QDP::QDP_isInitialized())
    QDP::QDP_initialize(&argc, &argv);
#endif
  
  Program_Parameters params;

  try {
    params = parse_input_arguments(argc, argv);
  } catch (std::exception &err) {
    QDPIO::cerr << "Error: " << err.what() << "\n" << std::endl;
    return 1;
  }

  // Set lattice dimensions
  /*
   * This needs to be determined from the scidac file.
   * Also, the other parameters in the scidac file
   * should be compared against the currently used QDP++
   */
  QDP::multi1d<int> nrow(Nd);
  for(int i=0; i<Nd; i++)
    nrow[i]=4;
  QDP::Layout::setLattSize(nrow);
  QDP::Layout::create();

  // following chroma file /mainprogs/main/purgaug.cc
  // Start up the config
  multi1d<LatticeColorMatrix> u(Nd);
  {
    QDP::XMLReader file_xml;
    QDP::XMLReader config_xml;
    std::string cfg_file = params.input_file;
    Chroma::readGauge(file_xml, config_xml, u, cfg_file, QDP::QDPIO_SERIAL);
  }


  QDPIO::cout << "QDP config imported, average plaquette:     "
            << fastsum::average_plaquette(u)
            << std::endl;

  try {
    mma::dump_lattice(params.output_file, u);
  } catch (std::exception &err) {
    QDPIO::cerr << err.what() << std::endl;
    return 1;
  }

}

