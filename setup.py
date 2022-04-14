from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'RAIChU'
LONG_DESCRIPTION = 'An interactive PKS/NRPS visualizer.'

setup(
    name="raichu",
    version=VERSION,
    author="Barbara Terlouw",
    author_email="barbara.r.terlouw@gmail.com",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    package_data={"": ["*.png", "*.svg"]},
    install_requires=[],
    scripts=["bin/raichu"],)