
/*
 * Created: 05-06-2017
 * Modified: Mon 05 Jun 2017 17:27:51 BST
 * Author: Jonas R. Glesaaen (jonas@glesaaen.com)
 */

#ifndef QDP_STREAM_REDIRECT_HPP
#define QDP_STREAM_REDIRECT_HPP

#include <qdp_stdio.h>
#include <iostream>
#include <fstream>

namespace fastsum {

class QDP_Stream_Redirect
{
public:
  QDP_Stream_Redirect(std::string new_cout_filename)
      : QDP_Stream_Redirect{new_cout_filename, new_cout_filename}
  {
  }

  QDP_Stream_Redirect(std::string new_cout_filename,
                      std::string new_cerr_filename)
      : cout_ofs{new_cout_filename}, cerr_ofs{new_cerr_filename}
  {
    if (!cout_ofs)
      throw std::runtime_error{"Cannot open \"" + new_cout_filename +
                               "\" for writing"};

    if (!cerr_ofs)
      throw std::runtime_error{"Cannot open \"" + new_cerr_filename +
                               "\" for writing"};

    QDP::QDPIO::cout.init(&cout_ofs);
    QDP::QDPIO::cerr.init(&cerr_ofs);
  }

  ~QDP_Stream_Redirect()
  {
    QDP::QDPIO::cout.init(&std::cout);
    QDP::QDPIO::cerr.init(&std::cerr);
  }

private:
  std::ofstream cout_ofs;
  std::ofstream cerr_ofs;
};

} // namespace fastsum

#endif /* QDP_STREAM_REDIRECT_HPP */
