
/*
 * Created: 06-06-2017
 * Modified: Wed 07 Jun 2017 16:46:31 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */

#include <openqcd_archive.hpp>
#include <stdexcept>

extern "C" {
#include "archive.h"
}

namespace fastsum {

/* Exports the current configuration stored in udfld() to a file named
 * "filename". See the OpenQCD documentation to see which parameters this
 * includes.
 */
void export_openqcd_config(std::string filename)
{
  if (!check_bc(0.0) or chs_ubnd(1))
    throw std::runtime_error {"Unexpeced boundary values for OpenQCD config"};

  // export_cnfg does not change the string, but it is not declared const
  export_cnfg(const_cast<char*>(filename.c_str()));
}

} // namespace fastsum 
