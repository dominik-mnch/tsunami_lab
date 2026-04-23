{ pkgs ? import <nixpkgs> {} }:

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
  ];
}
