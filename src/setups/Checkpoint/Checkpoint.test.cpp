/**
 * @author Magdalena Schwarzkopf
 * @author Dominik Münch
 *
 * @section DESCRIPTION
 * Unit tests for the Checkpoint setup.
 **/
#include <catch2/catch.hpp>
#include "Checkpoint.h"
#include "../../io/NetCdf.h"
#include <cstdio>
#include <filesystem>
#include <netcdf.h>
#include <string>
#include <vector>

/**
 * Creates a minimal checkpoint.nc with simulation parameters as global attributes.
 */
static void createCheckpointFile( std::string const & i_path,
                                  int i_nx, int i_ny,
                                  float i_dx, float i_dy,
                                  float i_originX, float i_originY,
                                  float i_endTime,
                                  float i_simTime = 0.0f,
                                  int i_k = 1,
                                  int i_solverMode = 1,
                                  std::string const & i_propagation = "2d",
                                  std::string const & i_setup = "artificial_tsunami_2d" ) {
  int l_ncId = -1;
  REQUIRE( nc_create( i_path.c_str(), NC_NETCDF4 | NC_CLOBBER, &l_ncId ) == NC_NOERR );
  REQUIRE( nc_put_att_int(   l_ncId, NC_GLOBAL, "nx",       NC_INT,   1, &i_nx      ) == NC_NOERR );
  REQUIRE( nc_put_att_int(   l_ncId, NC_GLOBAL, "ny",       NC_INT,   1, &i_ny      ) == NC_NOERR );
  REQUIRE( nc_put_att_int(   l_ncId, NC_GLOBAL, "k",        NC_INT,   1, &i_k       ) == NC_NOERR );
  REQUIRE( nc_put_att_int(   l_ncId, NC_GLOBAL, "solver_mode", NC_INT, 1, &i_solverMode ) == NC_NOERR );
  REQUIRE( nc_put_att_float( l_ncId, NC_GLOBAL, "dx",       NC_FLOAT, 1, &i_dx      ) == NC_NOERR );
  REQUIRE( nc_put_att_float( l_ncId, NC_GLOBAL, "dy",       NC_FLOAT, 1, &i_dy      ) == NC_NOERR );
  REQUIRE( nc_put_att_float( l_ncId, NC_GLOBAL, "origin_x", NC_FLOAT, 1, &i_originX ) == NC_NOERR );
  REQUIRE( nc_put_att_float( l_ncId, NC_GLOBAL, "origin_y", NC_FLOAT, 1, &i_originY ) == NC_NOERR );
  REQUIRE( nc_put_att_float( l_ncId, NC_GLOBAL, "end_time", NC_FLOAT, 1, &i_endTime ) == NC_NOERR );
  REQUIRE( nc_put_att_float( l_ncId, NC_GLOBAL, "sim_time", NC_FLOAT, 1, &i_simTime ) == NC_NOERR );
  REQUIRE( nc_put_att_text( l_ncId, NC_GLOBAL, "propagation", i_propagation.size(), i_propagation.c_str() ) == NC_NOERR );
  REQUIRE( nc_put_att_text( l_ncId, NC_GLOBAL, "setup", i_setup.size(), i_setup.c_str() ) == NC_NOERR );
  REQUIRE( nc_close( l_ncId ) == NC_NOERR );
}

/**
 * Creates a minimal solution.nc with two time steps.
 * The last fully committed step (t=5.0) holds h/hu/hv/b values we can verify.
 * A partial third step has height written but time is fill-value (simulates a crash).
 */
static void createSolutionFile( std::string const & i_path,
                                tsunami_lab::t_idx i_nx,
                                tsunami_lab::t_idx i_ny ) {
  int l_ncId = -1;
  REQUIRE( nc_create( i_path.c_str(), NC_NETCDF4 | NC_CLOBBER, &l_ncId ) == NC_NOERR );

  int l_dimTimeId = -1, l_dimYId = -1, l_dimXId = -1;
  REQUIRE( nc_def_dim( l_ncId, "time", NC_UNLIMITED,      &l_dimTimeId ) == NC_NOERR );
  REQUIRE( nc_def_dim( l_ncId, "y",    (std::size_t)i_ny, &l_dimYId    ) == NC_NOERR );
  REQUIRE( nc_def_dim( l_ncId, "x",    (std::size_t)i_nx, &l_dimXId    ) == NC_NOERR );

  int l_dimsYX[2]  = { l_dimYId, l_dimXId };
  int l_dimsTYX[3] = { l_dimTimeId, l_dimYId, l_dimXId };

  int l_varTimeId = -1, l_varHId = -1, l_varHuId = -1, l_varHvId = -1, l_varBId = -1;
  int l_varXId = -1, l_varYId = -1;
  REQUIRE( nc_def_var( l_ncId, "time",       NC_FLOAT, 1, &l_dimTimeId, &l_varTimeId ) == NC_NOERR );
  REQUIRE( nc_def_var( l_ncId, "x",          NC_FLOAT, 1, &l_dimXId,    &l_varXId    ) == NC_NOERR );
  REQUIRE( nc_def_var( l_ncId, "y",          NC_FLOAT, 1, &l_dimYId,    &l_varYId    ) == NC_NOERR );
  REQUIRE( nc_def_var( l_ncId, "height",     NC_FLOAT, 3, l_dimsTYX,    &l_varHId    ) == NC_NOERR );
  REQUIRE( nc_def_var( l_ncId, "momentum_x", NC_FLOAT, 3, l_dimsTYX,    &l_varHuId   ) == NC_NOERR );
  REQUIRE( nc_def_var( l_ncId, "momentum_y", NC_FLOAT, 3, l_dimsTYX,    &l_varHvId   ) == NC_NOERR );
  REQUIRE( nc_def_var( l_ncId, "bathymetry", NC_FLOAT, 2, l_dimsYX,     &l_varBId    ) == NC_NOERR );
  REQUIRE( nc_enddef( l_ncId ) == NC_NOERR );

  // ── time step 0 (t = 2.5) ──────────────────────────────────────────────────
  {
    float l_t = 2.5f;
    std::size_t l_s0[1] = {0}, l_c1[1] = {1};
    REQUIRE( nc_put_vara_float( l_ncId, l_varTimeId, l_s0, l_c1, &l_t ) == NC_NOERR );

    // h: all 1.0
    std::vector<float> l_h( i_nx * i_ny, 1.0f );
    std::size_t l_s3[3] = {0,0,0}, l_c3[3] = {1,(std::size_t)i_ny,(std::size_t)i_nx};
    REQUIRE( nc_put_vara_float( l_ncId, l_varHId,  l_s3, l_c3, l_h.data() ) == NC_NOERR );
    std::vector<float> l_hu( i_nx * i_ny, 0.5f );
    REQUIRE( nc_put_vara_float( l_ncId, l_varHuId, l_s3, l_c3, l_hu.data() ) == NC_NOERR );
    std::vector<float> l_hv( i_nx * i_ny, 0.25f );
    REQUIRE( nc_put_vara_float( l_ncId, l_varHvId, l_s3, l_c3, l_hv.data() ) == NC_NOERR );
  }

  // ── time step 1 (t = 5.0) — this is the committed checkpoint ──────────────
  {
    // data first, time last (commit marker protocol)
    std::size_t l_s3[3] = {1,0,0}, l_c3[3] = {1,(std::size_t)i_ny,(std::size_t)i_nx};
    // h: cell (iy,ix) = iy*nx + ix + 10
    std::vector<float> l_h( i_nx * i_ny );
    for( tsunami_lab::t_idx l_iy = 0; l_iy < i_ny; l_iy++ )
      for( tsunami_lab::t_idx l_ix = 0; l_ix < i_nx; l_ix++ )
        l_h[l_iy * i_nx + l_ix] = static_cast<float>( l_iy * i_nx + l_ix + 10 );
    REQUIRE( nc_put_vara_float( l_ncId, l_varHId, l_s3, l_c3, l_h.data() ) == NC_NOERR );

    std::vector<float> l_hu( i_nx * i_ny, 3.0f );
    REQUIRE( nc_put_vara_float( l_ncId, l_varHuId, l_s3, l_c3, l_hu.data() ) == NC_NOERR );

    std::vector<float> l_hv( i_nx * i_ny, 4.0f );
    REQUIRE( nc_put_vara_float( l_ncId, l_varHvId, l_s3, l_c3, l_hv.data() ) == NC_NOERR );

    float l_t = 5.0f;
    std::size_t l_s1[1] = {1}, l_c1[1] = {1};
    REQUIRE( nc_put_vara_float( l_ncId, l_varTimeId, l_s1, l_c1, &l_t ) == NC_NOERR );
  }

  // ── time step 2 — partial write simulating crash (height written, time missing) ──
  {
    std::size_t l_s3[3] = {2,0,0}, l_c3[3] = {1,(std::size_t)i_ny,(std::size_t)i_nx};
    std::vector<float> l_h( i_nx * i_ny, 999.0f );
    REQUIRE( nc_put_vara_float( l_ncId, l_varHId, l_s3, l_c3, l_h.data() ) == NC_NOERR );
    // time NOT written for step 2 — simulates crash before commit
  }

  // ── static bathymetry ──────────────────────────────────────────────────────
  {
    std::vector<float> l_b( i_nx * i_ny );
    for( tsunami_lab::t_idx l_i = 0; l_i < i_nx * i_ny; l_i++ )
      l_b[l_i] = static_cast<float>( l_i ) * -10.0f;
    REQUIRE( nc_put_var_float( l_ncId, l_varBId, l_b.data() ) == NC_NOERR );
  }

  REQUIRE( nc_close( l_ncId ) == NC_NOERR );
}

TEST_CASE( "Checkpoint setup loads simulation parameters from checkpoint.nc.", "[Checkpoint]" ) {
  std::string l_dir    = "test_cp_params_dir";
  std::string l_cpFile = l_dir + "/checkpoint.nc";
  std::string l_solFile = l_dir + "/solution.nc";
  std::filesystem::create_directories( l_dir );
  std::remove( l_cpFile.c_str() );
  std::remove( l_solFile.c_str() );

  constexpr tsunami_lab::t_idx l_nx = 3;
  constexpr tsunami_lab::t_idx l_ny = 2;

  createCheckpointFile( l_cpFile, l_nx, l_ny, 500.0f, 500.0f, 100.0f, 200.0f, 5.0f );
  createSolutionFile(   l_solFile, l_nx, l_ny );

  tsunami_lab::setups::Checkpoint l_cp( l_cpFile );

  REQUIRE( l_cp.getNx()      == l_nx    );
  REQUIRE( l_cp.getNy()      == l_ny    );
  REQUIRE( l_cp.getDx()      == Approx( 500.0 ) );
  REQUIRE( l_cp.getDy()      == Approx( 500.0 ) );
  REQUIRE( l_cp.getOriginX() == Approx( 100.0 ) );
  REQUIRE( l_cp.getOriginY() == Approx( 200.0 ) );
  REQUIRE( l_cp.getEndTime() == Approx( 5.0   ) );
  REQUIRE( l_cp.getK() == 1 );
  REQUIRE( l_cp.getSolverMode() == 1 );
  REQUIRE( l_cp.getPropagation() == "2d" );
  REQUIRE( l_cp.getSetup() == "artificial_tsunami_2d" );

  std::remove( l_cpFile.c_str() );
  std::remove( l_solFile.c_str() );
  std::filesystem::remove( l_dir );
}

TEST_CASE( "Checkpoint setup expands coarsened solution data using k.", "[Checkpoint]" ) {
  constexpr tsunami_lab::t_idx l_nx = 4;
  constexpr tsunami_lab::t_idx l_ny = 2;
  constexpr int l_k = 2;
  std::string l_dir = "test_cp_coarse_dir";
  std::string l_cpFile = l_dir + "/checkpoint.nc";
  std::string l_solFile = l_dir + "/solution.nc";
  std::filesystem::create_directories( l_dir );
  std::remove( l_cpFile.c_str() );
  std::remove( l_solFile.c_str() );

  createCheckpointFile( l_cpFile,
                        l_nx,
                        l_ny,
                        10.0f,
                        20.0f,
                        0.0f,
                        0.0f,
                        1.0f,
                        1.0f,
                        l_k,
                        0,
                        "2d",
                        "circular_dam_break_2d" );

  int l_ncId = -1;
  REQUIRE( nc_create( l_solFile.c_str(), NC_NETCDF4 | NC_CLOBBER, &l_ncId ) == NC_NOERR );
  int l_dimTimeId = -1, l_dimYId = -1, l_dimXId = -1;
  REQUIRE( nc_def_dim( l_ncId, "time", NC_UNLIMITED, &l_dimTimeId ) == NC_NOERR );
  REQUIRE( nc_def_dim( l_ncId, "y", l_ny / l_k, &l_dimYId ) == NC_NOERR );
  REQUIRE( nc_def_dim( l_ncId, "x", l_nx / l_k, &l_dimXId ) == NC_NOERR );
  int l_dimsYX[2] = { l_dimYId, l_dimXId };
  int l_dimsTYX[3] = { l_dimTimeId, l_dimYId, l_dimXId };
  int l_varTimeId = -1, l_varHId = -1, l_varHuId = -1, l_varHvId = -1, l_varBId = -1;
  REQUIRE( nc_def_var( l_ncId, "time", NC_FLOAT, 1, &l_dimTimeId, &l_varTimeId ) == NC_NOERR );
  REQUIRE( nc_def_var( l_ncId, "height", NC_FLOAT, 3, l_dimsTYX, &l_varHId ) == NC_NOERR );
  REQUIRE( nc_def_var( l_ncId, "momentum_x", NC_FLOAT, 3, l_dimsTYX, &l_varHuId ) == NC_NOERR );
  REQUIRE( nc_def_var( l_ncId, "momentum_y", NC_FLOAT, 3, l_dimsTYX, &l_varHvId ) == NC_NOERR );
  REQUIRE( nc_def_var( l_ncId, "bathymetry", NC_FLOAT, 2, l_dimsYX, &l_varBId ) == NC_NOERR );
  REQUIRE( nc_enddef( l_ncId ) == NC_NOERR );

  std::size_t l_s3[3] = { 0, 0, 0 };
  std::size_t l_c3[3] = { 1, l_ny / l_k, l_nx / l_k };
  std::vector<float> l_h = { 10.0f, 20.0f };
  std::vector<float> l_hu = { 1.0f, 2.0f };
  std::vector<float> l_hv = { 3.0f, 4.0f };
  std::vector<float> l_b = { -10.0f, -20.0f };
  REQUIRE( nc_put_vara_float( l_ncId, l_varHId, l_s3, l_c3, l_h.data() ) == NC_NOERR );
  REQUIRE( nc_put_vara_float( l_ncId, l_varHuId, l_s3, l_c3, l_hu.data() ) == NC_NOERR );
  REQUIRE( nc_put_vara_float( l_ncId, l_varHvId, l_s3, l_c3, l_hv.data() ) == NC_NOERR );
  std::size_t l_s2[2] = { 0, 0 };
  std::size_t l_c2[2] = { l_ny / l_k, l_nx / l_k };
  REQUIRE( nc_put_vara_float( l_ncId, l_varBId, l_s2, l_c2, l_b.data() ) == NC_NOERR );
  float l_t = 1.0f;
  std::size_t l_s1[1] = { 0 };
  std::size_t l_c1[1] = { 1 };
  REQUIRE( nc_put_vara_float( l_ncId, l_varTimeId, l_s1, l_c1, &l_t ) == NC_NOERR );
  REQUIRE( nc_close( l_ncId ) == NC_NOERR );

  tsunami_lab::setups::Checkpoint l_cp( l_cpFile );
  REQUIRE( l_cp.getK() == l_k );
  REQUIRE( l_cp.getSolverMode() == 0 );
  REQUIRE( l_cp.getPropagation() == "2d" );
  REQUIRE( l_cp.getSetup() == "circular_dam_break_2d" );

  for( tsunami_lab::t_idx l_iy = 0; l_iy < l_ny; l_iy++ ) {
    for( tsunami_lab::t_idx l_ix = 0; l_ix < l_nx; l_ix++ ) {
      tsunami_lab::t_real l_x = ( l_ix + 0.5f ) * 10.0f;
      tsunami_lab::t_real l_y = ( l_iy + 0.5f ) * 20.0f;
      tsunami_lab::t_real l_expectedBase = l_ix < 2 ? 10.0f : 20.0f;
      REQUIRE( l_cp.getHeight( l_x, l_y ) == Approx( l_expectedBase ) );
      REQUIRE( l_cp.getMomentumX( l_x, l_y ) == Approx( l_expectedBase / 10.0f ) );
      REQUIRE( l_cp.getMomentumY( l_x, l_y ) == Approx( l_expectedBase / 10.0f + 2.0f ) );
      REQUIRE( l_cp.getBathymetry( l_x, l_y ) == Approx( -l_expectedBase ) );
    }
  }

  std::remove( l_cpFile.c_str() );
  std::remove( l_solFile.c_str() );
  std::filesystem::remove( l_dir );
}

TEST_CASE( "Checkpoint setup loads last committed state and ignores partial step.", "[Checkpoint]" ) {
  constexpr tsunami_lab::t_idx l_nx = 3;
  constexpr tsunami_lab::t_idx l_ny = 2;
  constexpr float l_dx = 1000.0f;
  constexpr float l_dy = 1000.0f;
  constexpr float l_originX = 0.0f;
  constexpr float l_originY = 0.0f;

  // The Checkpoint constructor derives the solution path as sibling "solution.nc".
  // For the test we write to a temp dir so both files share the same directory.
  std::string l_dir      = "test_cp_dir";
  std::filesystem::create_directories( l_dir );
  std::string l_cpPath   = l_dir + "/checkpoint.nc";
  std::string l_solPath  = l_dir + "/solution.nc";
  std::remove( l_cpPath.c_str() );
  std::remove( l_solPath.c_str() );

  createCheckpointFile( l_cpPath, l_nx, l_ny, l_dx, l_dy, l_originX, l_originY, 5.0f );

  // Write solution.nc directly into the temp dir
  {
    int l_ncId = -1;
    REQUIRE( nc_create( l_solPath.c_str(), NC_NETCDF4 | NC_CLOBBER, &l_ncId ) == NC_NOERR );

    int l_dimTimeId = -1, l_dimYId = -1, l_dimXId = -1;
    REQUIRE( nc_def_dim( l_ncId, "time", NC_UNLIMITED,  &l_dimTimeId ) == NC_NOERR );
    REQUIRE( nc_def_dim( l_ncId, "y",    l_ny,           &l_dimYId    ) == NC_NOERR );
    REQUIRE( nc_def_dim( l_ncId, "x",    l_nx,           &l_dimXId    ) == NC_NOERR );

    int l_dimsYX[2]  = { l_dimYId, l_dimXId };
    int l_dimsTYX[3] = { l_dimTimeId, l_dimYId, l_dimXId };

    int l_varTimeId = -1, l_varHId = -1, l_varHuId = -1, l_varHvId = -1, l_varBId = -1;
    int l_varXId = -1, l_varYId = -1;
    REQUIRE( nc_def_var( l_ncId, "time",       NC_FLOAT, 1, &l_dimTimeId, &l_varTimeId ) == NC_NOERR );
    REQUIRE( nc_def_var( l_ncId, "x",          NC_FLOAT, 1, &l_dimXId,    &l_varXId    ) == NC_NOERR );
    REQUIRE( nc_def_var( l_ncId, "y",          NC_FLOAT, 1, &l_dimYId,    &l_varYId    ) == NC_NOERR );
    REQUIRE( nc_def_var( l_ncId, "height",     NC_FLOAT, 3, l_dimsTYX,    &l_varHId    ) == NC_NOERR );
    REQUIRE( nc_def_var( l_ncId, "momentum_x", NC_FLOAT, 3, l_dimsTYX,    &l_varHuId   ) == NC_NOERR );
    REQUIRE( nc_def_var( l_ncId, "momentum_y", NC_FLOAT, 3, l_dimsTYX,    &l_varHvId   ) == NC_NOERR );
    REQUIRE( nc_def_var( l_ncId, "bathymetry", NC_FLOAT, 2, l_dimsYX,     &l_varBId    ) == NC_NOERR );
    REQUIRE( nc_enddef( l_ncId ) == NC_NOERR );

    // step 0 — t=2.5, h=1, hu=0.5, hv=0.25
    {
      std::size_t l_s3[3] = {0,0,0}, l_c3[3] = {1,l_ny,l_nx};
      std::vector<float> l_hBuf( l_nx*l_ny, 1.0f );
      REQUIRE( nc_put_vara_float( l_ncId, l_varHId,  l_s3, l_c3, l_hBuf.data() ) == NC_NOERR );
      std::vector<float> l_huBuf( l_nx*l_ny, 0.5f );
      REQUIRE( nc_put_vara_float( l_ncId, l_varHuId, l_s3, l_c3, l_huBuf.data() ) == NC_NOERR );
      std::vector<float> l_hvBuf( l_nx*l_ny, 0.25f );
      REQUIRE( nc_put_vara_float( l_ncId, l_varHvId, l_s3, l_c3, l_hvBuf.data() ) == NC_NOERR );
      float l_t = 2.5f;
      std::size_t l_s1[1]={0}, l_c1[1]={1};
      REQUIRE( nc_put_vara_float( l_ncId, l_varTimeId, l_s1, l_c1, &l_t ) == NC_NOERR );
    }

    // step 1 (committed) — data first, time last
    {
      std::size_t l_s3[3] = {1,0,0}, l_c3[3] = {1,l_ny,l_nx};
      // h[iy][ix] = iy*nx + ix + 10
      std::vector<float> l_hBuf( l_nx*l_ny );
      for( tsunami_lab::t_idx l_iy = 0; l_iy < l_ny; l_iy++ )
        for( tsunami_lab::t_idx l_ix = 0; l_ix < l_nx; l_ix++ )
          l_hBuf[l_iy*l_nx+l_ix] = static_cast<float>( l_iy*l_nx+l_ix+10 );
      REQUIRE( nc_put_vara_float( l_ncId, l_varHId,  l_s3, l_c3, l_hBuf.data() ) == NC_NOERR );
      std::vector<float> l_huBuf( l_nx*l_ny, 3.0f );
      REQUIRE( nc_put_vara_float( l_ncId, l_varHuId, l_s3, l_c3, l_huBuf.data() ) == NC_NOERR );
      std::vector<float> l_hvBuf( l_nx*l_ny, 4.0f );
      REQUIRE( nc_put_vara_float( l_ncId, l_varHvId, l_s3, l_c3, l_hvBuf.data() ) == NC_NOERR );
      float l_t = 5.0f;
      std::size_t l_s1[1]={1}, l_c1[1]={1};
      REQUIRE( nc_put_vara_float( l_ncId, l_varTimeId, l_s1, l_c1, &l_t ) == NC_NOERR );
    }

    // step 2 — partial (crash sim): height written, time NOT written
    {
      std::size_t l_s3[3] = {2,0,0}, l_c3[3] = {1,l_ny,l_nx};
      std::vector<float> l_hBuf( l_nx*l_ny, 999.0f );
      REQUIRE( nc_put_vara_float( l_ncId, l_varHId, l_s3, l_c3, l_hBuf.data() ) == NC_NOERR );
    }

    // bathymetry: b[i] = i * -10
    {
      std::vector<float> l_bBuf( l_nx*l_ny );
      for( tsunami_lab::t_idx l_i = 0; l_i < l_nx*l_ny; l_i++ )
        l_bBuf[l_i] = static_cast<float>( l_i ) * -10.0f;
      REQUIRE( nc_put_var_float( l_ncId, l_varBId, l_bBuf.data() ) == NC_NOERR );
    }

    REQUIRE( nc_close( l_ncId ) == NC_NOERR );
  }

  tsunami_lab::setups::Checkpoint l_cp( l_cpPath );

  // ── simulation parameters ──────────────────────────────────────────────────
  REQUIRE( l_cp.getNx()      == l_nx );
  REQUIRE( l_cp.getNy()      == l_ny );
  REQUIRE( l_cp.getDx()      == Approx( l_dx ) );
  REQUIRE( l_cp.getDy()      == Approx( l_dy ) );
  REQUIRE( l_cp.getOriginX() == Approx( l_originX ) );
  REQUIRE( l_cp.getOriginY() == Approx( l_originY ) );
  REQUIRE( l_cp.getEndTime() == Approx( 5.0 ) );

  // ── state loaded from step 1 (not the partial step 2) ─────────────────────
  // Cell centres: x = originX + (ix + 0.5) * dx,  y = originY + (iy + 0.5) * dy
  for( tsunami_lab::t_idx l_iy = 0; l_iy < l_ny; l_iy++ ) {
    for( tsunami_lab::t_idx l_ix = 0; l_ix < l_nx; l_ix++ ) {
      tsunami_lab::t_real l_x = l_originX + ( l_ix + 0.5f ) * l_dx;
      tsunami_lab::t_real l_y = l_originY + ( l_iy + 0.5f ) * l_dy;

      float l_expectedH = static_cast<float>( l_iy * l_nx + l_ix + 10 );
      float l_expectedB = static_cast<float>( l_iy * l_nx + l_ix ) * -10.0f;

      REQUIRE( l_cp.getHeight(     l_x, l_y ) == Approx( l_expectedH ) );
      REQUIRE( l_cp.getMomentumX(  l_x, l_y ) == Approx( 3.0 ) );
      REQUIRE( l_cp.getMomentumY(  l_x, l_y ) == Approx( 4.0 ) );
      REQUIRE( l_cp.getBathymetry( l_x, l_y ) == Approx( l_expectedB ) );
    }
  }

  std::remove( l_cpPath.c_str() );
  std::remove( l_solPath.c_str() );
  std::filesystem::remove( l_dir );
}
