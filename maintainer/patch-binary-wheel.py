"""Patch binary wheel to allow executable binaries.

Following the patch in https://github.com/pypa/auditwheel/pull/443/
"""

import base64
import hashlib
import logging
import shutil
import sys
import tempfile
from pathlib import Path

logging.basicConfig()
log = logging.getLogger("binary-patch")
log.setLevel(logging.INFO)

WRAPPER = """
#!python
import os
import sys
import sysconfig
if __name__ == "__main__":
    os.execv(
        os.path.join(sysconfig.get_path("platlib"), {binary_path!r}),
        sys.argv,
    )
"""


def hash(path: Path):
    """Compute file hash."""
    return (
        base64.urlsafe_b64encode(hashlib.sha256(path.read_bytes()).digest())
        .decode()
        .rstrip("=")
    )


def record(tmpd):
    """Generate a RECORD file.

    Rehash the content of the wheel.
    """
    lines = []
    for path in Path(tmpd).glob("**/*"):
        if not path.is_file():
            continue

        path_ = path.relative_to(tmpd)
        if path.name == "RECORD":
            lines.append(f"{path_},,")
            continue

        hash_ = hash(path)
        size = path.stat().st_size
        lines.append(f"{path_},sha256={hash_},{size}")

    return "\n".join(lines)


def patch(wheel: Path):
    """Patch wheel to add wrapper script."""
    with tempfile.TemporaryDirectory() as tmpd:
        shutil.unpack_archive(wheel, tmpd, format="zip")
        scripts = Path(tmpd) / "pineappl_cli.scripts"
        scripts.mkdir()
        installed = next(Path(tmpd).glob("*.data")) / "scripts" / "pineappl"
        shutil.move(installed, scripts)
        binary = scripts / "pineappl"
        installed.write_text(WRAPPER.format(binary_path=str(binary.relative_to(tmpd))))
        next(Path(tmpd).glob("**/RECORD")).write_text(record(tmpd))
        shutil.make_archive(wheel.name, format="zip", base_dir=".", root_dir=tmpd)
        shutil.move(wheel.name + ".zip", wheel)
        log.info(f"Patched {wheel} in-place")


def main():
    """Patch wheels passed as arguments."""
    wheels = [Path(p) for p in sys.argv[1:]]
    log.info("Processing the following wheels:\n" + "\n".join(str(w) for w in wheels))
    for wheel in wheels:
        patch(wheel)


if __name__ == "__main__":
    main()
