/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Checkpoint setup: delegates all netCDF I/O to NetCdf::readCheckpoint.
 **/
#include "Checkpoint.h"
#include "../../io/NetCdf.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

tsunami_lab::setups::Checkpoint::Checkpoint( std::string const & i_checkpointFile ) {
  tsunami_lab::io::NetCdf::CheckpointData l_data =
    tsunami_lab::io::NetCdf::readCheckpoint( i_checkpointFile );

  m_nx      = l_data.nx;
  m_ny      = l_data.ny;
  m_dx      = l_data.dx;
  m_dy      = l_data.dy;
  m_originX = l_data.originX;
  m_originY = l_data.originY;
  m_simTime = l_data.simTime;
  m_endTime = l_data.endTime;
  m_h       = std::move( l_data.h  );
  m_hu      = std::move( l_data.hu );
  m_hv      = std::move( l_data.hv );
  m_b       = std::move( l_data.b  );
}

tsunami_lab::t_idx tsunami_lab::setups::Checkpoint::getIndex( t_real i_x, t_real i_y ) const {
  t_idx l_ix = static_cast<t_idx>( std::max( t_real(0),
                 std::min( t_real(m_nx - 1), std::floor( (i_x - m_originX) / m_dx ) ) ) );
  t_idx l_iy = static_cast<t_idx>( std::max( t_real(0),
                 std::min( t_real(m_ny - 1), std::floor( (i_y - m_originY) / m_dy ) ) ) );
  return l_iy * m_nx + l_ix;
}

tsunami_lab::t_real tsunami_lab::setups::Checkpoint::getHeight( t_real i_x, t_real i_y ) const {
  return m_h[ getIndex( i_x, i_y ) ];
}

tsunami_lab::t_real tsunami_lab::setups::Checkpoint::getBathymetry( t_real i_x, t_real i_y ) const {
  return m_b[ getIndex( i_x, i_y ) ];
}

tsunami_lab::t_real tsunami_lab::setups::Checkpoint::getMomentumX( t_real i_x, t_real i_y ) const {
  return m_hu[ getIndex( i_x, i_y ) ];
}

tsunami_lab::t_real tsunami_lab::setups::Checkpoint::getMomentumY( t_real i_x, t_real i_y ) const {
  return m_hv[ getIndex( i_x, i_y ) ];
}

