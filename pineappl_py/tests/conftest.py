import subprocess

import pytest


class PDF:
    def xfxQ2(self, pid, x, q2):
        if pid in range(-6, 6):
            return x * (1 - x)
        else:
            return 0.0

    def xfxQ(self, pid, x, q):
        return self.xfxQ2(pid, x, q**2)

    def alphasQ(self, q):
        return 1.0

    # Define the Toy Polarized PDF set
    def polarized_pdf(self, pid, x, q2):
        return 2.0

    # Define the Toy Unpolarized PDF set
    def unpolarized_pdf(self, pid, x, q2):
        return 1.0


@pytest.fixture
def pdf():
    return PDF()


@pytest.fixture
def download_objects(tmp_path_factory):
    def _download_fk(objname: str) -> None:
        download_dir = tmp_path_factory.mktemp("data")
        file_path = download_dir / f"{objname}"
        args = [
            "wget",
            "--no-verbose",
            "--no-clobber",
            "-P",
            f"{download_dir}",
            f"https://data.nnpdf.science/pineappl/test-data/{objname}",
        ]

        try:
            _ = subprocess.run(
                args,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return file_path
        except OSError as error:
            msg = f"Failed to execute the command {args}."
            raise EnvironmentError(msg) from error

    return _download_fk
