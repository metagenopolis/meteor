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
