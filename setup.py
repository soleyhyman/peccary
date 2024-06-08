#!/usr/bin/env python3
#
# This file is part of the peccary technique package
#
#
#
#
import os
import setuptools

# find here
here = os.path.abspath(os.path.dirname(__file__))

# ---- Perform setup                                                        ----
# get long description
with open("README.md", "r") as fh:
    long_description = fh.read()

# setuptools.setup(use_scm_version=True)
setuptools.setup(
    name="peccary",
    version="0.1",
    author="Sóley Hyman",
    author_email="soleyhyman@arizona.edu",
    description="Permutation Entropy and Statistical Complexity Analysis for Astrophysics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(include=['peccary','peccary.*']),
    python_requires='>=3.6',
    install_requires=["numpy >= 1.14","scipy >= 0.19","matplotlib"],
    package_data={"": ["README.md","LICENSE"]},
    license='GPL-3.0',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Software Development",
        "Topic :: Software Development :: Libraries :: Python Modules"]
)