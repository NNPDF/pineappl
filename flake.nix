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
      # Common arguments can be set here to avoid repeating them later
      commonArgs = {
        inherit src;
        strictDeps = true;
      };

      cargoArtifacts = craneLib.buildDepsOnly commonArgs;
      individualCrateArgs =
        commonArgs
        // {
          inherit cargoArtifacts;
          inherit
            ((builtins.fromTOML (builtins.readFile (src
                    + "/Cargo.toml")))
              .workspace
              .package)
            version
            ;
        };

      fileSetForCrate = crate:
        lib.fileset.toSource {
          root = ./.;
          fileset = lib.fileset.unions [
            ./Cargo.toml
            ./Cargo.lock
            crate
          ];
        };

      pineappl = craneLib.buildPackage (individualCrateArgs
        // {
          pname = "pineappl";
          cargoExtraArgs = "-p pineappl";
          src = fileSetForCrate ./pineappl;
        });
      cli = craneLib.buildPackage (individualCrateArgs
        // {
          pname = "pineappl_cli";
          cargoExtraArgs = "-p pineappl_cli";
          src = fileSetForCrate ./pineappl_cli;
          buildInputs = with pkgs; [pkg-config lhapdf];
        });
      # TODO: build with maturin
      py = null;
    in
      {
        packages = {
          default = cli;
          cli = cli;
          lib = pineappl;
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
