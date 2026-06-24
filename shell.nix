{ pkgs ? import <nixpkgs> { config = { allowUnfree = true; }; } }:

pkgs.mkShell {
  buildInputs = with pkgs; [
    python3
    python3Packages.pandas
    python3Packages.matplotlib
    python3Packages.imageio
    python3Packages.numpy
    python3Packages.pillow
    python3Packages.sphinx
    python3Packages.sphinx-rtd-theme
    doxygen
    scons
    gmt
    netcdf
    paraview
    cudaPackages.cudatoolkit
    cudaPackages.cuda_nvcc
  ];

  shellHook = ''
    export CUDA_PATH=${pkgs.cudaPackages.cudatoolkit}
    export PATH=$CUDA_PATH/bin:$PATH
    export CPATH=$CUDA_PATH/include:$CPATH
    export LD_LIBRARY_PATH=$CUDA_PATH/lib64:$CUDA_PATH/lib:$LD_LIBRARY_PATH
    echo "CUDA toolkit loaded"
  '';
}
