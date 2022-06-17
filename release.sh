NEW_VERSION=$1

echo "Releasing pineappl $NEW_VERSION"

# only update main crate and bindings
# not plugins
sed -E -i 's/(^version = ").*"/\1'$NEW_VERSION'"/' \
  pineappl/Cargo.toml \
  pineappl_capi/Cargo.toml \
  pineappl_cli/Cargo.toml \
  pineappl_py/Cargo.toml

# always replace pineappl version wherever it is a dependency
find . -name Cargo.toml -exec \
  sed -E -i 's/(pineappl .* version = ").*"/\1'$NEW_VERSION'"/' {} +

echo "Versions updated, check git ;)"
