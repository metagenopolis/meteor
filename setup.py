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

"""Setup meteor"""
# pylint: disable=deprecated-module
# pylint: disable=redefined-builtin
# import os
from setuptools import setup, find_packages
# from setuptools.command.install import install
# from distutils.command.build import build
# from subprocess import call

# METEOR_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "meteor", "src")

# class MeteorBuild(build):
#     """Compile for meteor binaries"""
#     def run(self):
#         # run original build code
#         build.run(self)
#         # build meteor
#         build_path = os.path.abspath(self.build_temp)
#         cmd = [
#             "make",
#             "OUT=" + build_path
#         ]

#         target_files = [os.path.join(build_path, "meteor-counter"),
#                         os.path.join(build_path, "meteor-profiler")]

#         def compile():
#             call(cmd, cwd=METEOR_PATH)

#         self.execute(compile, [], "Compiling meteor")

#         # copy resulting tool to library build folder
#         self.mkpath(self.build_scripts)

#         if not self.dry_run:
#             for target in target_files:
#                 self.copy_file(target, self.build_scripts)


# class MeteorInstall(install):
#     """Install process for meteor binaries"""
#     def initialize_options(self):
#         install.initialize_options(self)
#         self.build_scripts = None

#     def finalize_options(self):
#         install.finalize_options(self)
#         self.set_undefined_options("build", ("build_scripts", "build_scripts"))

#     def run(self):
#         # run original install code
#         install.run(self)
#         # install Meteor executables
#         self.copy_tree(self.build_scripts, self.install_scripts)



with open("README.md", encoding="UTF-8") as f:
    readme = f.read()
setup(name="meteor",
      version="3.3",
      license="GPLv3",
      description="A plateform for quantitative metagenomic profiling of complex ecosystems.",
      long_description=readme,
      author="Amine Ghozlane",
      author_email="amine.ghozlane@pasteur.fr",
      platforms= ["Linux", "Unix", "Darwin", "Windows"],
      install_requires=["pysam", "pandas"],
      packages=find_packages(),
      classifiers = [
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    package_data = {
        "meteor": ["*.ini"]
    },
    entry_points={
        "console_scripts": [
           "meteor=meteor.meteor:main",
        ]
    },
    # cmdclass={
    #     "build": MeteorBuild,
    #     "install": MeteorInstall,
    # }
)
