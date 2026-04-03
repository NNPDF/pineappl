# Maintainer-specific tools

- `Containerfile`: description to generate the container
- `download-cli-dependencies.sh`: Bash script to install dependencies needed by
  the CLI and used when generating a container and when making releases
- `build-container.sh`: Bash script to build the container
- `download-test-data.sh`: downloads test data used for testing and code
  coverage
- `generate-coverage.sh`: generates code coverage locally
- `make-release.sh`: this script is supposed to be run by the maintainer to
  create a new release. The only argument is the version number
- `reinterpolate-eko.py`: Python script to reinterpolate an EKO to used with
  PineAPPL's `evolve` subcommand
- `test-cli-export.sh`: downloads PineAPPL grids from [Ploughshare] and tests
  PineAPPL's CLI `export` subcommand with them
- `test-cli-import.sh`: downloads APPLgrids and fastNLO tables from
  [Ploughshare] and [Mitov's ttbar page] and tests PineAPPL's CLI `import`
  subcommand with them

[Ploughshare]: https://ploughshare.web.cern.ch/ploughshare/
[Mitov's ttbar page]: https://www.precision.hep.phy.cam.ac.uk/results/ttbar-fastnlo/
