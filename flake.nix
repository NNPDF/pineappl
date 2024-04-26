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
      # The following could be done automatically by Crane, but it will look for a
      # top-level package name as well (not present in PineAPPL). Thus, doing it
      # manually will save a Crane warning.
      version =
        (builtins.fromTOML (builtins.readFile (src
              + "/Cargo.toml")))
        .workspace
        .package
        .version;
      commonArgs = {
        inherit src version;
        strictDeps = false;
        buildInputs = with pkgs; [lhapdf];
        nativeBuildInputs = with pkgs; [pkg-config];
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

      defaultName = "pineappl";
      cargoArtifacts = craneLib.buildDepsOnly (commonArgs // {pname = defaultName;});
      cli = craneLib.buildPackage (
        commonArgs
        // {
          pname = defaultName;
          cargoExtraArgs = "-p pineappl_cli";
          inherit cargoArtifacts;
          src = fileSetForCrate ./pineappl_cli;
        }
      );
      # TODO: build with maturin
      py =
        (craneLib.buildPackage (commonArgs
          // {
            pname = "pineappl-py";
            inherit cargoArtifacts;
          }))
        .overrideAttrs (old: {
          nativeBuildInputs = old.nativeBuildInputs ++ [pkgs.maturin];
          buildPhase =
            old.buildPhase
            + ''
              maturin build --offline --target-dir ./target
            '';
          installPhase =
            old.installPhase
            + ''
              cp target/wheels/pineappl $out/
            '';
        });
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
