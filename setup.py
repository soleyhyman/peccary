import os
from setuptools import setup, find_packages

# find here
here = os.path.abspath(os.path.dirname(__file__))

# get long description
with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name='peccary',
    version='0.1.0',
    description='Package for identifying regular, complex, and stochastic behavior in timeseries',
    long_description=long_description,
    long_description_content_type="text/rst",
    url='https://github.com/soleyhyman/peccary',
    author='Sóley Hyman',
    author_email='soleyhyman@arizona.edu',
    packages=find_packages(include=["peccary","peccary.*"]),
    python_requires='>=3.8',
    install_requires=["numpy >= 1.14","scipy >= 0.19","matplotlib"],
    package_data={"": ["README.rst","LICENSE"]},
    license='GPL-3.0',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)