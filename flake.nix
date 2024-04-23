{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.05";
    flake-utils.url = "github:numtide/flake-utils";
    devenv.url = "github:cachix/devenv";
    crane = {
      url = "github:ipetkov/crane";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    fenix.url = "github:nix-community/fenix";
    fenix.inputs.nixpkgs.follows = "nixpkgs";
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
    crane,
    devenv,
    ...
  } @ inputs:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = nixpkgs.legacyPackages.${system};
    in {
      packages = let
        craneLib = crane.lib.${system};
        pineappl = craneLib.buildPackage {
          src = craneLib.cleanCargoSource (craneLib.path ./.);
          strictDeps = true;
          # buildInputs = [ # Add additional build inputs here ]
        };
      in {
        default = {
          packages.default = pineappl;
          packages.pineappl-py = null;
        };
      };

      devShells.default = let
        pwd = builtins.getEnv "PWD";
        prefix = "${pwd}/target/prefix";
        lhapath = "${prefix}/share/LHAPDF";
      in
        devenv.lib.mkShell {
          inherit inputs pkgs;
          modules = [
            ({config, ...}: {
              packages = with pkgs; [
                maturin
                (lhapdf.override {
                  python =
                    config.languages.python.package;
                })
              ];

              env = {
                PREFIX = prefix;
                LHAPDF_DATA_PATH = lhapath;
              };
              enterShell = ''
                # update path before entering the shell, when Nix packages updates
                # already happened
                export PATH=${prefix}/bin:$PATH
                mkdir -p ${lhapath}
              '';

              languages.python = {
                enable = true;
                venv.enable = true;
              };
              languages.rust = {
                enable = true;
                channel = "stable";
              };
            })
          ];
        };
    });

  nixConfig = {
    extra-trusted-public-keys = "devenv.cachix.org-1:w1cLUi8dv3hnoSPGAuibQv+f9TZLr6cv/Hm9XgU50cw=";
    extra-substituters = "https://devenv.cachix.org";
  };
}
