{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
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
    in let
      inherit (pkgs) lib;
      craneLib = crane.lib.${system};

      src = craneLib.cleanCargoSource (craneLib.path ./.);
      commonArgs = {
        inherit src;
        strictDeps = false;
        buildInputs = with pkgs; [lhapdf];
        nativeBuildInputs = with pkgs; [pkg-config];
      };

      individualCrateArgs =
        commonArgs
        // {
          inherit
            ((builtins.fromTOML (builtins.readFile (src
                    + "/Cargo.toml")))
              .workspace
              .package)
            version
            ;
          doCheck = false;
        };

      fileSetForCrate = crate:
        lib.fileset.toSource {
          root = ./.;
          fileset = lib.fileset.unions [
            ./Cargo.toml
            ./Cargo.lock
            crate
            # TODO: avoid passing all the folders
            # (if strictly needed, at least avoid doing it explicitly)
            ./pineappl
            ./pineappl_applgrid
            ./pineappl_capi
            ./pineappl_cli
            ./pineappl_fastnlo
            ./pineappl_py
            ./xtask
          ];
        };

      cli = craneLib.buildPackage (individualCrateArgs
        // {
          pname = "pineappl";
          cargoExtraArgs = "-p pineappl_cli";
          src = fileSetForCrate ./pineappl_cli;
        });
      # TODO: build with maturin
      py = null;
    in
      {
        packages = {
          default = cli;
          cli = cli;
          pineappl-py = py;
        };

        apps.default = cli;
      }
      // {
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
