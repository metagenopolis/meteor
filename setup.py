import os
from setuptools import setup, find_packages
from setuptools.command.install import install
from distutils.command.build import build
from subprocess import call

METEOR_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'meteor', 'src')

class MeteorBuild(build):
    def run(self):
        # run original build code
        build.run(self)
        # build meteor
        build_path = os.path.abspath(self.build_temp)
        cmd = [
            'make',
            'OUT=' + build_path
        ]

        target_files = [os.path.join(build_path, 'meteor-counter'),
                        os.path.join(build_path, 'meteor-profiler')]

        def compile():
            call(cmd, cwd=METEOR_PATH)

        self.execute(compile, [], 'Compiling meteor')

        # copy resulting tool to library build folder
        self.mkpath(self.build_scripts)

        if not self.dry_run:
            for target in target_files:
                self.copy_file(target, self.build_scripts)


class MeteorInstall(install):
    def initialize_options(self):
        install.initialize_options(self)
        self.build_scripts = None

    def finalize_options(self):
        install.finalize_options(self)
        self.set_undefined_options('build', ('build_scripts', 'build_scripts'))

    def run(self):
        # run original install code
        install.run(self)
        # # install Meteor executables
        # target_files = [os.path.join(self.build_lib, 'meteor-counter'),
        #                 os.path.join(self.build_lib, 'meteor-profiler')]
        # for target in target_files:
        #         self.copy_file(target, self.install_scripts)
        self.copy_tree(self.build_scripts, self.install_scripts)
        # self.copy_tree(self.build_lib, self.install_scripts)



with open('README.md') as f:
    readme = f.read()
setup(name='meteor',
      version='3.3',
      license='GPLv3',
      description='A plateform for quantitative metagenomic profiling of complex ecosystems.',
      long_description=readme,
      author='Amine Ghozlane',
      author_email='amine.ghozlane@pasteur.fr',
      platforms= ['Linux', 'Unix', 'Darwin', 'Windows'],
      install_requires=[],
      packages=find_packages(),
      # package_dir = {'meteor': 'meteor'},
      classifiers = [
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    package_data = {
        'meteor': ['*.ini']
    },
    entry_points={
        'console_scripts': [
           'meteor=meteor.meteor:main',
        ]
    },
    cmdclass={
        'build': MeteorBuild,
        'install': MeteorInstall,
    }
)
