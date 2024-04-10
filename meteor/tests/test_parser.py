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

"""Test parser main object"""

# pylint: disable=redefined-outer-name
from ..session import Component
from ..parser import Parser
from pathlib import Path
import pytest


@pytest.fixture
def parser_standard(datadir: Path) -> Parser:
    meteor = Component
    meteor.profile_dir = datadir / "profiles"
    module_file = datadir / "module_def.tsv"
    return Parser(module_file=module_file)


def test_load_modules(parser_standard: Parser) -> None:
    module_dict = parser_standard.load_modules(parser_standard.module_file)
    expected_output = {
        "MF0005": "K01639 ((K00885+K01788),K13967)",
        "MF0006": "((K01941+K01457),(K01428+K01429+K01430),K14048)",
        "MF0009": "K01667",
        "MF0011": "(K00812,K00813,K11358)",
        "MF0013": "(K00260,K15371,K00261,K00262) K00109 ((K01039+K01040)) K01726 K01615",
        "MF0015": "K01846 K04835",
        "MF0016": "(K00281,(K00282+K00283)) K00605 K00382",
        "MF0018": "((K00318+K00294),K13821)",
        "MF0019": "K01777 ((K10793+K10794+K10795+K10796))",
        "MF0026": "(K01697,K10150) K01758",
        "MF0027": "(K01758,K01760)",
        "MF0029": "((K00060+K00639),(K01620+K00132),(K01620+K04072),(K01620+K04073))",
        "MF0033": "((K00812+K01011),(K00813+K01011),(K11358+K01011))",
        "MF0034": "((K00265+K00266))",
        "MF0035": "(K01583,K01584,K01585,K02626) ((K10536+K12251))",
        "MF0041": "K01745 K01712 K01468 (K01479,(K05603+K01458),(K00603+K13990))",
    }
    assert module_dict == expected_output


def test_clean_module(parser_standard: Parser) -> None:
    mod_def = parser_standard.clean_module("K1 K2")
    assert mod_def == "K1+K2"
    mod_def = parser_standard.clean_module("(K1 K2)")
    assert mod_def == "K1+K2"
    mod_def = parser_standard.clean_module("K1 (K2,K3)")
    assert mod_def == "(K1+(K2,K3))"
    mod_def = parser_standard.clean_module("K1-K7 (K2,K3)")
    assert mod_def == "(K1+(K2,K3))"
    mod_def = parser_standard.clean_module("K1 (K2,K3+K4)+K5")
    assert mod_def == "(K1+(K2,K3+K4)+K5)"


def test_replace_submod(parser_standard: Parser) -> None:
    mod_def = "K01,K02 M03+K04"
    mod_dict = {"M01": "K01", "M03": "K0123+K0124"}
    mod_def = parser_standard.replace_submod(mod_def, mod_dict)
    assert mod_def == "K01,K02 K0123+K0124+K04"


def test_find_all_alt(parser_standard: Parser) -> None:
    mod_dict = {"M01": "K01", "M03": "K0123+K0124"}
    real_alt = parser_standard.find_all_alt("K01 K02+K03+(K04,K05)", mod_dict)
    true_alt = [{"K01", "K02", "K03", "K04"}, {"K01", "K02", "K03", "K05"}]
    assert all(x in true_alt for x in real_alt)
    assert len(real_alt) == len(true_alt)
    real_alt = parser_standard.find_all_alt("(K01,K02) K03+(K04,K05)", mod_dict)
    true_alt = [
        {"K01", "K03", "K04"},
        {"K01", "K03", "K05"},
        {"K02", "K03", "K04"},
        {"K02", "K03", "K05"},
    ]
    assert all(x in true_alt for x in real_alt)
    assert len(real_alt) == len(true_alt)
    real_alt = parser_standard.find_all_alt(
        "((K16154+K16155),(K16157+K16158),K08684)", mod_dict
    )
    true_alt = [
        {"K16154", "K16155"},
        {"K16157", "K16158"},
        {"K08684"},
    ]
    assert all(x in true_alt for x in real_alt)
    assert len(real_alt) == len(true_alt)
    real_alt = parser_standard.find_all_alt(
        "(K13811,(K00957+K00956)) ((K00394+K00395)) ((K11180+K11181))", mod_dict
    )
    true_alt = [
        {"K13811", "K00394", "K00395", "K11180", "K11181"},
        {"K00957", "K00956", "K00394", "K00395", "K11180", "K11181"},
    ]
    assert all(x in true_alt for x in real_alt)
    assert len(real_alt) == len(true_alt)


def test_execute(parser_standard: Parser) -> None:
    parser_standard.execute()
    expected_output = {
        "MF0005": [{"K01639", "K00885", "K01788"}, {"K01639", "K13967"}],
        "MF0006": [{"K01941", "K01457"}, {"K01428", "K01429", "K01430"}, {"K14048"}],
        "MF0009": [{"K01667"}],
        "MF0011": [{"K00812"}, {"K00813"}, {"K11358"}],
        "MF0013": [
            {"K00260", "K00109", "K01039", "K01040", "K01726", "K01615"},
            {"K15371", "K00109", "K01039", "K01040", "K01726", "K01615"},
            {"K00261", "K00109", "K01039", "K01040", "K01726", "K01615"},
            {"K00262", "K00109", "K01039", "K01040", "K01726", "K01615"},
        ],
        "MF0015": [{"K01846", "K04835"}],
        "MF0016": [
            {"K00281", "K00605", "K00382"},
            {"K00282", "K00283", "K00605", "K00382"},
        ],
        "MF0018": [{"K00318", "K00294"}, {"K13821"}],
        "MF0019": [{"K01777", "K10793", "K10794", "K10795", "K10796"}],
        "MF0026": [{"K01697", "K01758"}, {"K10150", "K01758"}],
        "MF0027": [{"K01758"}, {"K01760"}],
        "MF0029": [
            {"K00060", "K00639"},
            {"K01620", "K00132"},
            {"K01620", "K04072"},
            {"K01620", "K04073"},
        ],
        "MF0033": [{"K00812", "K01011"}, {"K00813", "K01011"}, {"K11358", "K01011"}],
        "MF0034": [{"K00265", "K00266"}],
        "MF0035": [
            {"K01583", "K10536", "K12251"},
            {"K01584", "K10536", "K12251"},
            {"K01585", "K10536", "K12251"},
            {"K02626", "K10536", "K12251"},
        ],
        "MF0041": [
            {"K01745", "K01712", "K01468", "K01479"},
            {"K01745", "K01712", "K01468", "K05603", "K01458"},
            {"K01745", "K01712", "K01468", "K00603", "K13990"},
        ],
    }
    assert len(expected_output) == len(parser_standard.module_dict_alt)
    assert all(x in expected_output for x in parser_standard.module_dict_alt)
    assert all(x in parser_standard.module_dict_alt for x in expected_output)
