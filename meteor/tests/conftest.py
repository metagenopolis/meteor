import pytest
import shutil
from pathlib import Path


@pytest.fixture(autouse=True)
def cleanup(tmp_path: Path, request: pytest.FixtureRequest) -> None:
    # This fixture will run before and after each test
    def teardown():
        # Remove the temporary directory and its contents
        for item in tmp_path.iterdir():
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()

    request.addfinalizer(teardown)


FAKE_FREEBAYES = r'''
import gzip
import os
import sys

args = sys.argv[1:]
if args == ["--version"]:
    print("version:  v1.3.6")
    sys.exit(0)
if os.environ.get("FAKE_FREEBAYES_FAIL"):
    sys.stderr.write("fake freebayes failure\n")
    sys.exit(1)
source = os.environ.get("FAKE_FREEBAYES_VCF")
if source:
    with gzip.open(source, "rt") as vcf:
        sys.stdout.write(vcf.read())
else:
    reference = args[args.index("-f") + 1]
    print("##fileformat=VCFv4.2")
    print('##INFO=<ID=RO,Number=1,Type=Integer,Description="ref count">')
    print('##INFO=<ID=AO,Number=A,Type=Integer,Description="alt count">')
    with open(reference + ".fai") as fai:
        for line in fai:
            name, length = line.split("\t")[:2]
            print(f"##contig=<ID={name},length={length}>")
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
'''


@pytest.fixture
def fake_freebayes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> pytest.MonkeyPatch:
    """A `freebayes` executable on PATH: prints FAKE_FREEBAYES_VCF (gzipped VCF)
    when set, otherwise a VCF header without record; fails when
    FAKE_FREEBAYES_FAIL is set. Lets the variant calling run without freebayes."""
    import os
    import sys

    bindir = tmp_path / "fake_bin"
    bindir.mkdir()
    script = bindir / "freebayes"
    script.write_text(f"#!{sys.executable}\n{FAKE_FREEBAYES}")
    script.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.delenv("FAKE_FREEBAYES_VCF", raising=False)
    monkeypatch.delenv("FAKE_FREEBAYES_FAIL", raising=False)
    return monkeypatch
