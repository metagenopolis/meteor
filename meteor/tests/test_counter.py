# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Test counter builder main objects"""

# pylint: disable=redefined-outer-name
from ..session import Component
from ..counter import Counter
from pathlib import Path
import pytest


@pytest.fixture
def counter_best(tmp_path: Path):
    meteor = Component
    meteor.ref_dir = tmp_path
    meteor.ref_name = "test"
    meteor.threads = 1
    return Counter(meteor, counting_type="best", mapping_type="end-to-end", trim=80,
                   alignment_number=10000, counting_only=False, mapping_only=False)