"""
To run package in development mode, run `pip install -e .` in this folder
"""

from setuptools import find_packages, setup

setup(
    name="vectorseq",
    packages=find_packages(),
    version="0.3.0",
    description="VECTORseq Data Analysis",
    author="Victoria Cheung, Philip Chung",
    license="MIT",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3 :: Only",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
    ],
)
