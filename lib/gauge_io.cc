/*
 * Routines associated with Chroma gauge config IO
 * Copied from chroma project lib/io/gaugeio.cc
 */

#include "gauge_io.h"

namespace Chroma {


// Read a Chroma propagator
/*
 * \param file_xml     xml reader holding config info ( Modify )
 * \param record_xml   xml reader holding config info ( Modify )
 * \param u            gauge configuration ( Modify )
 * \param file         path ( Read )
 * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    
void readGauge(XMLReader& file_xml,
	       XMLReader& record_xml, 
	       multi1d<LatticeColorMatrix>& u, 
	       const std::string& file, 
	       QDP_serialparallel_t serpar)
{

  /* 
   * This constructor assumes that the lattice
   * dimensions have already been defined.
   */
  QDPFileReader to(file_xml,file,serpar);

  /* 
   * This is problematic - the size should
   * come from the read - a resize. Currently, QDPIO does not 
   * support this
   */
  read(to,record_xml,u);
  if (to.bad())
  {
    QDPIO::cerr << __func__ << ": error reading file " << file << std::endl;
    QDP_abort(1);
  }
 
  close(to);
}


// Write a Gauge field in QIO format
/*
 * \param file_xml    xml reader holding config info ( Modify )
 * \param record_xml  xml reader holding config info ( Modify )
 * \param u           gauge configuration ( Modify )
 * \param file        path ( Read )
 * \param volfmt      either QDPIO_SINGLEFILE, QDPIO_MULTIFILE ( Read )
 * \param serpar      either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    
void writeGauge(XMLBufferWriter& file_xml,
		XMLBufferWriter& record_xml, 
		const multi1d<LatticeColorMatrix>& u, 
		const std::string& file, 
		QDP_volfmt_t volfmt, 
		QDP_serialparallel_t serpar)
{
  QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
  if (to.bad())
  {
    QDPIO::cerr << __func__ << ": error writing file " << file << std::endl;
    QDP_abort(1);
  }

  write(to, record_xml, u);         // Write in native precision
  close(to);
}


}  // end namespace Chroma
