import os
# from sys import platform
from setuptools import setup
from setuptools.command.install import install
from distutils.command.build import build
from subprocess import call

METEOR_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'meteor', 'src')
print(METEOR_PATH)
class MeteorBuild(build):
    def run(self):
        # run original build code
        build.run(self)
        # build meteor
        build_path = os.path.abspath(self.build_temp)
        print(build_path)
        cmd = [
            'make',
            'OUT=' + build_path
        ]

        # targets = ['python']
        # cmd.extend(targets)

        # if platform == 'darwin':
        #     target_path = 'OSX64_PYTHON'
        # else:
        #     target_path = 'UNIX_PYTHON'

        target_files = [os.path.join(build_path, 'meteor-counter'), 
                        os.path.join(build_path, 'meteor-profiler')]

        def compile():
            call(cmd, cwd=METEOR_PATH)

        self.execute(compile, [], 'Compiling meteor')

        # copy resulting tool to library build folder
        self.mkpath(self.build_lib)

        if not self.dry_run:
            for target in target_files:
                self.copy_file(target, self.build_lib)


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

        # install Meteor executables
        self.copy_tree(self.build_lib, self.install_lib)


with open('README.rst') as f:
    readme = f.read()
setup(name='meteor',
      version='4.3',
      license='GPLv3',
      description='A plateform for quantitative metagenomic profiling of complex ecosystems.',
      long_description=readme,
      author='Amine Ghozlane',
      author_email='amine.ghozlane@pasteur.fr',
      platforms= ['Linux', 'Unix', 'Darwin', 'Windows'],
      package_dir={'meteor':'meteor'},
      classifiers = [
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
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