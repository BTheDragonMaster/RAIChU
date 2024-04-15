from setuptools import setup, find_packages


VERSION = '1.0.2'
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
                      'biopython==1.76',
                      'pikachu-chem>=1.1.1',
                      'timeout-decorator==0.5.0'],
    scripts=["bin/raichu"],)
