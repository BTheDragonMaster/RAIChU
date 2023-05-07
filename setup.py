from setuptools import setup, find_packages

VERSION = '0.0.2'
DESCRIPTION = 'RAIChU'
LONG_DESCRIPTION = 'A package to automatically draw natural product biosynthesis pathways.'

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
                      'pikachu-chem>=1.0.13',
                      'paras'],
    scripts=["bin/raichu"],)
