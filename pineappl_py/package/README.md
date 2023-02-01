# Packaging python interface

In order to compile wheels to distribute some requirements have to be met:

- `linux`: the compilation process has to be run in a
  [`manylinux`](https://github.com/pypa/manylinux) compliant environment, for
  this reason a suitable container image is provided (see
  [published packages](https://github.com/orgs/NNPDF/packages?repo_name=pineappl)
  and the respective [`Containerfile`](./Containerfile))

  Notice that the default container provided by
  [pypa](https://github.com/pypa/manylinux) is not sufficient, since it does not
  ship a C compiler (required to compile the `syn` crate).
- `macOS`: it just needs to be run in a macOS environment, see
  [publishing workflow](https://github.com/NNPDF/pineappl/tree/master/.github/workflows/wheels.yml)
- `windows`: it just needs to be run in a windows environment, see
  [publishing workflow](https://github.com/NNPDF/pineappl/tree/master/.github/workflows/wheels.yml)

## `maturin` container image

`maturin` image has its own version (following [semver](https://semver.org/)),
and:

- it is based on `manylinux2014_x86_64`
- build wheels for a range of CPython versions (the actual one depends on the
  `maturin` version inside the container)

### Using `maturin` to compile for `manylinux`

This is the easy part: you just need to download the
[image](https://github.com/NNPDF/pineappl/pkgs/container/maturin) and run with
your favorite container tool.

Here the explicit commands with `podman` [[1]](#docker)

```sh
podman pull ghcr.io/n3pdf/maturin:latest
podman run ghcr.io/n3pdf/maturin
podman cp <container-id>:root/pineappl/pineappl_py/target/wheels/ .
```

Now wheels are available outside the container and can be uploaded in your
favorite way.

#### Interactive use

If you want to use the container environment interactively, you need to provide
an alternative entry point:

```sh
podman run --entrypoint bash -it <your-image>
```

### Create a new `maturin` image

_Use case_: if a new rust or maturin version is released, it might be needed to
upgrade also those inside the `maturin` image (since they are pre-installed in
the image itself)

To upgrade the build instructions ([_Containerfile_](./Containerfile)):

- change `FROM` to choose a different manylinux base image (see the
  [official source](https://github.com/pypa/manylinux))
- change `ARG` for `MATURIN_TAG` to choose a different `maturin` version
- to change architecture both `FROM` and `ARG MATURIN_TAR=` have to be updated
  from `x86_64`
- `rust` version is always taken to be the latest one at each build (so
  rerunning the build **without changing** anything might generate a **different
  image**)

Once `Containerfile` has been updated then rerun:

```sh
# inside pineappl/pineappl_py/package
podman build -t ghcr.io/n3pdf/maturin .
podman tag <image-id> ghcr.io/n3pdf/maturin:<version>
# login to GitHub registry with user credentials (not organization), see [2]
echo ${PAT} | podman login ghcr.io -u <username> --password-stdin
# finally publish
podman push ghcr.io/n3pdf/maturin:<version>
# and publish the new latest (all layers already available, it's just an alias)
podman push ghcr.io/n3pdf/maturin:latest
```

<a name="docker">[1]</a>: In the following I will use `podman` as the container
runtime for the examples. To use `docker` instead, you can simply replace
`podman -> docker`, they have compatible subcommands
<a name="github-registry-docs">[2]</a>: official
[GitHub registry docs](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry)
