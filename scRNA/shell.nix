let
  pkgs = import <nixpkgs> {};

  buildRPackage = pkgs.rPackages.buildRPackage;

  rWrapper = pkgs.rWrapper.override {
    packages = with pkgs.rPackages; [
      ggplot2
      data_table
      MASS
      cowplot
      limma
    ];
  };
in
  pkgs.mkShell {
    buildInputs = [
      rWrapper
    ];
  }
