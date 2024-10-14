from setuptools import setup

setup(
    name="pycchem",
    version="1.0.0",
    description="A DFT post-processing package",
    maintainer="Mohan Shankar",
    maintainer_email="mjs7eek@virginia.edu",
    license="MIT",
    url="https://github.com/mohan-s1/pycchem",
    packages=[
        "pycchem",
        "pycchem.gaussian",
        "pycchem.vasp",
        "pycchem.misc_utils",
    ],
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "bs4",
        "pytest",
        "re",
        "scikit-learn",
    ],
    long_description="""Package
      includes utilties to calculate entropy from Gaussian and VASP frequency calculations.""",
)