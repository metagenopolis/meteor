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

"""Parse functional modules"""

from meteor.session import Session
from dataclasses import dataclass
import pandas as pd
from pathlib import Path
import re


@dataclass
class Parser(Session):
    """Parse functional modules."""

    module_file: Path

    def load_modules(self, module_file: Path) -> dict[str, str]:
        """Create the module dictionary starting from a tsv file

        :param module_file: path to the tab-separated module file (no header,
        four columns: id (e.g., MF0001), type (e.g., GMM), name (e.g., butyrate production),
        definition (e.g., K01+K02))
        """
        module = pd.read_table(module_file, names=["id", "type", "name", "definition"])
        # Remove trailing whitespaces if any
        module["definition"] = module["definition"].str.rstrip()
        module_dict = dict(zip(module.id, module.definition))
        return module_dict

    def clean_module(self, mod_def: str) -> str:
        """Clean a module definition so that there is only '+' or ',' or parenthesis.
        Components such as '-XXX' or '--' are removed.
        Whitespaces (e.g., different steps) are replaced by '+'.

        :param mod_def: a single module definition (e.g. K01+K02 (K03, K04))
        """
        # Remove ' -- '
        mod_def = re.sub(r"\s?--\s?", r" ", mod_def).strip()
        # Remove '-XXX'
        mod_def = re.sub(r"\-\w+", r"", mod_def)
        # Replace ' ' by '+'
        mod_def = mod_def.replace(" ", "+")
        # Add parenthesis around the whole mod_def (in case we have such case : KXXX,KYYY
        mod_def = "(" + mod_def + ")"
        # Remove single parenthesis (do it as long as the module_def is changed, ie, as long as the pattern is found
        # (the while and the flag are necessary for such patterns that could be found : ((KO1+KO2))
        flag = 1
        while flag:
            mod_def_new = re.sub(r"\(([\w+]+)\)", r"\1", mod_def)
            if mod_def_new == mod_def:
                flag = 0
            else:
                mod_def = mod_def_new
        return mod_def

    def replace_submod(self, mod_def: str, mod_dict: dict[str, str]) -> str:
        """Replace submodules from mod_def according to their definition in mod_dict

        :param mod_def: a module definition containing a submodule
        :param mod_dict: a module dictionary containing the definition of the submodule
        """
        mod_def = re.sub(
            r"\w+",
            lambda x: mod_dict[x.group()] if x.group() in mod_dict else x.group(),
            mod_def,
        )
        return mod_def

    def find_all_alt(self, mod_def: str, mod_dict: dict[str, str]) -> list[set[str]]:
        """Return the list of all alternatives for a module definition.
        Based on the following rules : '+' = complex, ',' = alternative.

        :param mod_def: a single module definition (e.g., "K01,K02+K03")
        :param mod_dict: a module dictionary for submodule solving
        """
        # Define an alternative
        alternative = re.compile(r"(.*)(\(([\w+]+,{1})+[\w+]+\))(.*)")
        list_alt = [mod_def]
        # Flag to know if some alternatives remain to be solved
        flag = True
        while flag:  # While alternatives remain
            flag = False  # At the beginning of each loop we suppose there are no more alternative
            new_list_alt = (
                []
            )  # Initialize a new list to keep alternatives that will be solved in the next loop
            for my_def in list_alt:
                my_def = self.clean_module(mod_def=my_def)
                my_def = self.replace_submod(mod_def=my_def, mod_dict=mod_dict)
                # Look for simple parenthesis (simple parenthesis = alternatives that should be solved)
                res = alternative.match(my_def)
                if res:
                    flag = (
                        True  # Alternatives were found so we suppose other may remain
                    )
                    # list of alternatives after resolving simple parenthesis
                    new_def = [
                        res.group(1) + k + res.group(4)
                        for k in res.group(2).strip("()").split(",")
                    ]
                    new_list_alt += new_def
                else:  # If the def has already been simplified at maximum
                    new_list_alt.append(my_def)
            list_alt = new_list_alt
        # Transform each alternative into set of KOs : module_dict_alt['M000x'] = [set(KO1, KO2), set(KO1, KO3), etc]
        return [set(k.split("+")) for k in set(list_alt)]

    def execute(self) -> None:
        "Parse a module definition file to get all the possible alternatives"
        # Load file
        module_dict = self.load_modules(self.module_file)
        # Get the list of alternatives for all modules
        self.module_dict_alt = {
            mod: self.find_all_alt(mod_def=my_def, mod_dict=module_dict)
            for (mod, my_def) in module_dict.items()
        }
