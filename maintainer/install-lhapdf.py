"""Install LHAPDF."""

import io
import logging
import os
import sys
import tarfile
import tempfile
from http.client import HTTPSConnection
from pathlib import Path
from subprocess import run

VERSION = "6.5.3"
TARBALL = f"LHAPDF-{VERSION}.tar.gz"
HOST = "lhapdf.hepforge.org"
PATH = f"/downloads/?f={TARBALL}"
SRC_DIR = f"LHAPDF-{VERSION}"


logging.basicConfig()
log = logging.getLogger("LHAPDF installer")
log.setLevel(logging.INFO)


def download(host: str, path: str) -> bytes:
    """Download content from remote server."""
    conn = HTTPSConnection(host)
    conn.request("GET", path)
    resp = conn.getresponse()
    return resp.read()


def install(path: Path):
    """Build from source and install."""

    def run_(cmd, *args, **kwargs):
        cmd = cmd.split() if isinstance(cmd, str) else cmd
        run(cmd, *args, cwd=path / SRC_DIR, **kwargs)

    run_("autoreconf -f -i")
    prefix = os.environ.get("PREFIX")
    prefix_ = ["--prefix", prefix] if prefix is not None else []
    run_(["./configure"] + prefix_, env={"PYTHON": sys.executable})
    log.info("Configured")
    run_("make clean")
    run_("make -j")
    log.info("Built")
    run_("make install")
    log.info("Installed")


def main():
    """Install LHAPDF."""
    tar = download(HOST, PATH)
    log.info(f"Downloaded {HOST}{PATH}")
    with tempfile.TemporaryDirectory() as tmpd:
        tarfile.open(fileobj=io.BytesIO(tar), mode="r:gz").extractall(tmpd)
        log.info(f"Extracted LHAPDF tarbal in {tmpd}")
        install(Path(tmpd))


if __name__ == "__main__":
    main()
