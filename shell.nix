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
    # On WSL the real NVIDIA driver (libcuda.so.1) lives in /usr/lib/wsl/lib.
    # It must precede the Nix CUDA libs so the runtime loads the host driver
    # instead of a stub, otherwise kernel launches fail with
    # "CUDA driver version is insufficient for CUDA runtime version".
    if [ -d /usr/lib/wsl/lib ]; then
      export LD_LIBRARY_PATH=/usr/lib/wsl/lib:$LD_LIBRARY_PATH
    fi
    echo "CUDA toolkit loaded"
  '';
}
