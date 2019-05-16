
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


/*
 * Read the scidac file name and the lattice size from
 * the purgaug xml file.
 * From Chroma file mainprogs/main/purgaug.cc
 * Stripped out the fields we don't need.
 */

 //! Holds params for Heat-bath
struct HBItrParams 
{
  QDP::multi1d<int> nrow;
  QDP::Real beta;
};
void read(QDP::XMLReader& xml_in, const std::string& path,
	  HBItrParams& hbitr_params,
	  std::string& cfg_file)
{
  try {
    QDP::XMLReader paramtop(xml_in, path);
    QDP::read(paramtop, "HBItr/nrow", hbitr_params.nrow);
    QDP::read(paramtop, "HBItr/GaugeAction/beta", hbitr_params.beta);
    QDP::read(paramtop, "Cfg/cfg_file", cfg_file);
  }
  catch(const std::string& e) {
    QDPIO::cerr << "Caught Exception reading HBControl: " << e << std::endl;
    QDP_abort(1);
  }
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

  /*
   * Input data file name and lattice dimensions
   * from purgaug output xml file.
   *
   * Ideally, one would get the dimensions directly from
   * the scidac (*.lime) file directly.  However,
   * the QIO library discover_dims_mode is triggered off
   * the lattice dimensions, which is a compile-time parameter
   * for QDP/QDP++.  Thus, the lattice dimensions should be checked
   * against the lime file.  
   * A better solution would be for discover_dims_mode to trigger
   * from lattice volume=0 or maybe one of the sizes=0.
   *
   * In addition, the QDP++ library is not set up to trigger
   * discover_dims_mode anyway.  Finally, it is not clear if 
   * discover_dims_mode has ever been tested.
   */
  HBItrParams hbitr_params;
  std::string cfg_file;
  try
    {
      QDP::XMLReader xml_in(params.input_file);
      read(xml_in, "/purgaug", hbitr_params, cfg_file);
    }
  catch( const std::string& e ) {
    QDPIO::cerr << "Caught Exception reading input XML: " << e << std::endl;
    QDP_abort(1);
  }
  catch( std::exception& e ) {
    QDPIO::cerr << "Caught standard library exception: " << e.what() << std::endl;
    QDP_abort(1);
  }
  catch(...) {
    QDPIO::cerr << "Caught unknown exception " << std::endl;
    QDP_abort(1);
  }

  // Set lattice dimensions
  /*
   * This could be determined from the scidac file.
   * Also, the other parameters in the scidac file
   * should be compared against the currently used QDP++
   */
  QDP::Layout::setLattSize(hbitr_params.nrow);
  QDP::Layout::create();

  // following chroma file /mainprogs/main/purgaug.cc
  // Start up the config
  multi1d<LatticeColorMatrix> u(Nd);
  {
    QDP::XMLReader file_xml;
    QDP::XMLReader config_xml;
    Chroma::readGauge(file_xml, config_xml, u, cfg_file, QDP::QDPIO_SERIAL);
  }


  QDPIO::cout << "QDP config imported, average plaquette:     "
            << fastsum::average_plaquette(u)
            << std::endl;

  try {
    mma::dump_lattice(params.output_file, u, hbitr_params.beta);
  } catch (std::exception &err) {
    QDPIO::cerr << err.what() << std::endl;
    return 1;
  }

}

