from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'RAIChU'
LONG_DESCRIPTION = 'An interactive PKS/NRPS visualizer.'

print(find_packages())

setup(
    name="raichu",
    version=VERSION,
    author="Barbara Terlouw",
    author_email="barbara.r.terlouw@gmail.com",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(exclude="build"),
    package_data={"": ["*.png", "*.svg", "*.txt"]},
    install_requires=['matplotlib',
                      'pikachu-chem',
                      'paras'],
    scripts=["bin/raichu"],)
