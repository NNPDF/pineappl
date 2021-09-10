# To run with podman
# $ podman pull quay.io/pypa/manylinux2014_x86_64
# $ podman run -it quay.io/pypa/manylinux2014_x86_64
# $ podman cp build.sh <CONTAINER-ID>:/
# # sh build.sh
#
# install c compiler
yum install gcc -y

# create a dedicated user
# useradd -m pineappl
# su - pineappl
cd /root

# podman cp pineappl sweet_gagarin:home/pineappl
git clone https://github.com/N3PDF/pineappl.git

# add prefix
mkdir -p local/bin
export PATH=${HOME}/local/bin:${HOME}/.cargo/bin:${PATH}

# install dependencies
# - rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
# - maturin
MATURIN_TAR='maturin-x86_64-unknown-linux-musl.tar.gz'
MATURIN_TAG='v0.11.3'
curl --remote-name -L \
  https://github.com/PyO3/maturin/releases/download/${MATURIN_TAG}/${MATURIN_TAR}
tar -xvzf ${MATURIN_TAR} --directory=local/bin/

cd pineappl/pineappl_py
maturin build --release #--strip
