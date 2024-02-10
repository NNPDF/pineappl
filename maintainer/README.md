# Maintainer-specific tools

- `pineappl-ci/Containerfile`: description to generate the container
- `pineappl-ci/script.sh`: Bash script to install all dependencies into the
  container
- `generate-coverage.sh`: generates code coverage locally
- `make-release.sh`: this script is supposed to be run by the maintainer to
  create a new release. The only argument is the version number
- `test-ploughshare-conversion.sh`: downloads APPLgrids and fastNLO tables from
  [Ploughshare] and tests PineAPPL's CLI `import` subcommand with them

[Ploughshare]: https://ploughshare.web.cern.ch/ploughshare/
