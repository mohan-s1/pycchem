from setuptools import setup, find_packages

setup(
    name="pycchem",
    version="1.0.0",
    description="A DFT post-processing package",
    maintainer="Mohan Shankar",
    maintainer_email="mjs7eek@virginia.edu",
    license="MIT",
    url="https://github.com/mohan-s1/pycchem",
    packages=find_packages(include=[
        "pycchem",
        "pycchem.*",
    ]),
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "bs4",
        "pytest",
        "scikit-learn",
    ],
    long_description="""Package
      includes utilties to calculate entropy from Gaussian and VASP frequency calculations.""",
)