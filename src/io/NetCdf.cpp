/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Interface for NetCDF.
 **/

#include <netcdf.h>
#include "NetCdf.h"

tsunami_lab::io::NetCdf::NetCdf ( t_real               i_dx,
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
									t_real       const * i_hv)
:
	m_dx(i_dx),
	m_dy(i_dy),
	m_originX(i_originX),
	m_originY(i_originY),
	m_nx(i_nx),
	m_ny(i_ny),
	m_stride(i_stride),
	m_time(i_time),
	m_h(i_h),
	m_b(i_b),
	m_hu(i_hu),
	m_hv(i_hv),
	
	// open netCDF file for writing, replace if it already exists
	m_file("solution.nc", netCDF::NcFile::replace){

    //the dimensions and vaiables will be defined here, and the header will be written to the file
	}

												

void tsunami_lab::io::NetCdf::writeTimeStep(t_real simTime) {

	// the actual writing of the data to the netCDF file will be implemented here, using the member variables for the data and the dimensions
    
}

