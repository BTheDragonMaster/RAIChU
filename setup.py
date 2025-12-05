from setuptools import setup, find_packages


VERSION = "1.0.4"
DESCRIPTION = "RAIChU"
LONG_DESCRIPTION = (
    "A package to automatically draw natural product biosynthesis pathways."
)

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
    python_requires=">=3.7, <3.10",
    install_requires=["matplotlib", "biopython==1.83", "pikachu-chem>=1.1.4", "timeout-decorator"],
)
