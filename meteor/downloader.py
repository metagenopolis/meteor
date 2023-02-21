from dataclasses import dataclass, field
from session import Session, Component
from referencebuilder import ReferenceBuilder
from configparser import ConfigParser
from pathlib import Path
from hashlib import md5
from tqdm import tqdm
from urllib.request import urlretrieve
import logging


"""
Download and index reference
"""

class TqdmUpTo(tqdm):
    """Provides `update_to(n)` which uses `tqdm.update(delta_n)`."""
    def update_to(self, b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            self.total = tsize
        return self.update(b * bsize - self.n)  # also sets self.n = b * bsize

@dataclass
class Downloader(Session):
    """Download and prepare catalogues"""
    # user_choice: list
    meteor: Component
    choice: str
    configuration_path: Path = field(default_factory=Path)
    catalogues_config: ConfigParser = field(default_factory=ConfigParser)

    def __post_init__(self)->None:
        self.configuration_path = Path(__file__).parent / "dataverse_inrae.ini"
        assert self.configuration_path.exists()
        self.catalogues_config.read_file(self.configuration_path.open("rt", encoding="UTF-8"))
        self.meteor.ref_dir.mkdir(exist_ok=True)

    def getmd5(self, catalog):
        logging.info("Compute md5sum")
        with open(catalog, "rb") as f:
            file_hash = md5()
            while chunk := f.read(8192):
                file_hash.update(chunk)
        return file_hash.hexdigest()

    def download_file(self, url, catalog_fasta):
        with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1,
              desc=url.split('/')[-1]) as t:  # all optional kwargs
            urlretrieve(url, filename=catalog_fasta,
                       reporthook=t.update_to, data=None)
            t.total = t.n

    def execute(self)->bool:
        try:
            # for choice in self.user_choice:
            logging.info(f"Download {self.choice} microbiome reference catalog")
            url = self.catalogues_config[self.choice]["nr_catalogue"]
            md5_expect = self.catalogues_config[self.choice]["md5"]
            catalog_fasta = self.meteor.ref_dir / self.catalogues_config[self.choice]["filename"]
            self.download_file(url, catalog_fasta)
            assert md5_expect == self.getmd5(catalog_fasta)
            reference_builder = ReferenceBuilder(self.meteor,
                catalog_fasta)
            reference_builder.execute()
        except AssertionError:
            logging.error(f"MD5sum of {catalog_fasta} has a different value than expected")
        return True
