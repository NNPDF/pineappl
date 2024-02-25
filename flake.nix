{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.05";
    systems.url = "github:nix-systems/default";
    devenv.url = "github:cachix/devenv";
    nixpkgs-python.url = "github:cachix/nixpkgs-python";
    fenix.url = "github:nix-community/fenix";
    fenix.inputs.nixpkgs.follows = "nixpkgs";
  };

  outputs = {
    self,
    nixpkgs,
    devenv,
    systems,
    ...
  } @ inputs: let
    forEachSystem = nixpkgs.lib.genAttrs (import systems);
  in {
    devShells =
      forEachSystem
      (system: let
        pkgs = nixpkgs.legacyPackages.${system};
        pwd = builtins.getEnv "PWD";
        prefix = "${pwd}/target/prefix";
        lhapath = "${prefix}/share/LHAPDF";
      in {
        default = devenv.lib.mkShell {
          inherit inputs pkgs;
          modules = [
            ({config, ...}: let
              path = config.system.path;
            in {
              packages = with pkgs; [
                maturin
                (lhapdf.override {
                  python =
                    config.languages.python.package;
                })
              ];

              env = {
                PREFIX = prefix;
                PATH = "${path}:${prefix}/bin";
                LHAPDF_DATA_PATH = lhapath;
              };
              enterShell = ''
                mkdir -p ${lhapath}
              '';

              languages.python = {
                enable = true;
                venv.enable = true;
                version = "3.11";
              };
              languages.rust = {
                enable = true;
                channel = "stable";
              };
            })
          ];
        };
      });
  };

  nixConfig = {
    extra-trusted-public-keys = "devenv.cachix.org-1:w1cLUi8dv3hnoSPGAuibQv+f9TZLr6cv/Hm9XgU50cw=";
    extra-substituters = "https://devenv.cachix.org";
  };
}
