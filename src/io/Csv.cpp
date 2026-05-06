/**
 * @author Alexander Breuer (alex.breuer AT uni-jena.de)
 *
 * @section DESCRIPTION
 * IO-routines for writing a snapshot as Comma Separated Values (CSV).
 **/
#include "Csv.h"

void tsunami_lab::io::Csv::write( t_real               i_dx,
                                  t_real               i_dy,
                                  t_real               i_originX,
                                  t_real               i_originY,
                                  t_idx                i_nx,
                                  t_idx                i_ny,
                                  t_idx                i_stride,
                                  t_real               i_time,
                                  t_real       const * i_h,
                                  t_real       const * i_b,
                                  t_real       const * i_hu,
                                  t_real       const * i_hv,
                                  std::ostream       & io_stream ) {
  // write the CSV header
  io_stream << "time,x,y";
  if( i_h  != nullptr ) io_stream << ",height";
  if( i_b  != nullptr ) io_stream << ",bathymetry";
  if( i_hu != nullptr ) io_stream << ",momentum_x";
  if( i_hv != nullptr ) io_stream << ",momentum_y";
  io_stream << "\n";

  // iterate over all cells
  for( t_idx l_iy = 0; l_iy < i_ny; l_iy++ ) {
    for( t_idx l_ix = 0; l_ix < i_nx; l_ix++ ) {
      // derive coordinates of cell center
      t_real l_posX = i_originX + (l_ix + 0.5) * i_dx;
      t_real l_posY = i_originY + (l_iy + 0.5) * i_dy;

      t_idx l_id = l_iy * i_stride + l_ix;

      // write data
      io_stream << i_time << "," << l_posX << "," << l_posY;
      if( i_h  != nullptr ) io_stream << "," << i_h[l_id];
      if( i_b  != nullptr ) io_stream << "," << i_b[l_id];
      if( i_hu != nullptr ) io_stream << "," << i_hu[l_id];
      if( i_hv != nullptr ) io_stream << "," << i_hv[l_id];
      io_stream << "\n";
    }
  }
  io_stream << std::flush;
}

bool tsunami_lab::io::Csv::readBathymetry( std::string const &          i_filePath,
                                           std::vector<t_real>        & o_x,
                                           std::vector<t_real>        & o_b ) {
  // clear output vectors
  o_x.clear();
  o_b.clear();

  // open file
  std::ifstream l_file( i_filePath );
  if( !l_file.is_open() ) {
    std::cerr << "Error: Could not open bathymetry file: " << i_filePath << std::endl;
    return false;
  }

  std::string l_line;

  // skip header line
  if( !std::getline( l_file, l_line ) ) {
    std::cerr << "Error: Bathymetry file is empty: " << i_filePath << std::endl;
    l_file.close();
    return false;
  }

  // read data lines
  while( std::getline( l_file, l_line ) ) {
    // skip empty lines
    if( l_line.empty() ) continue;

    std::istringstream l_stream( l_line );
    std::string l_token;
    std::vector<t_real> l_values;

    // parse comma-separated values
    while( std::getline( l_stream, l_token, ',' ) ) {
      try {
        l_values.push_back( std::stod( l_token ) );
      } catch ( ... ) {
        std::cerr << "Error: Could not parse value '" << l_token << "' as a number" << std::endl;
        l_file.close();
        return false;
      }
    }

    // check if we have the correct number of columns
    if( l_values.size() != 4 ) {
      std::cerr << "Error: Expected 4 columns (lat, long, x, height), but got " 
                << l_values.size() << std::endl;
      l_file.close();
      return false;
    }

    // store values in output vectors
    o_x.push_back( l_values[2] );
    o_b.push_back( l_values[3] );
  }

  l_file.close();
  return true;
}