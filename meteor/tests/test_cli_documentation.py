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

"""The README documents every option of every command"""

import argparse
import re
import sys
from pathlib import Path

import pytest

from .. import meteor as meteor_cli

README = Path(__file__).resolve().parents[2] / "README.md"


class _Parser(Exception):
    pass


@pytest.fixture
def subcommands(monkeypatch) -> dict[str, argparse.ArgumentParser]:
    """The subcommand parsers built by get_arguments"""

    def capture(self, *args, **kwargs):
        raise _Parser(self)

    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", capture)
    monkeypatch.setattr(sys, "argv", ["meteor", "test"])
    with pytest.raises(_Parser) as captured:
        meteor_cli.get_arguments()
    parser = captured.value.args[0]
    subparsers = next(a for a in parser._actions if isinstance(a, argparse._SubParsersAction))
    return {name: sub for name, sub in subparsers.choices.items() if name != "test"}


def _documented_options() -> dict[str, set[str]]:
    """{command: options} of the tables of the "Options of each command" section"""
    section = README.read_text(encoding="utf-8").split("## Options of each command", 1)[1]
    section = section.split("\n## ", 1)[0]
    documented = {}
    for block in re.split(r"^### meteor ", section, flags=re.M)[1:]:
        command, table = block.split("\n", 1)
        documented[command.strip()] = {
            flag
            for cell in re.findall(r"^\| `([^`]+)`", table, flags=re.M)
            for flag in re.findall(r"(?<![\w-])(--?[\w]+)", cell)
        }
    return documented


def test_readme_documents_every_option(subcommands) -> None:
    documented = _documented_options()
    assert set(documented) == set(subcommands)
    for command, sub in subcommands.items():
        options = {
            flag
            for action in sub._actions
            if not isinstance(action, argparse._HelpAction)
            for flag in action.option_strings
        }
        assert options == documented[command], command
