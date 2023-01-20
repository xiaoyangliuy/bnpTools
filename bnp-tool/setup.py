from setuptools import setup
from setuptools import find_packages
import os

this_dir = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(this_dir, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="gyl-tools",
    version="0.1",
    description="Data processing toolbox for BNP data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Yanqi Luo",
    author_email="yluo89@anl.gov",
    download_url="https://github.com/AdvancedPhotonSource/bnpTools/",
    license="MIT",
    install_requires=[
        "h5py",
        "matplotlib",
        "numpy",
        "ax-platform",
        "websockets",
        "tqdm",
        "numpy",
        "scikit-learn",
        "pystackreg"
    ],
    # extras_require={
    #     'model_saving': ['h5py'],
    #     'molecules': ['openbabel', 'rdkit'],
    #     'tensorflow': ['tensorflow>=2.1'],
    #     'tensorflow with gpu': ['tensorflow-gpu>=2.1'],
    # },
    packages=find_packages(),
    package_data={
        # "hardware": ["*.yaml", "*/*.yaml", "*/*/*.yaml", "*/*/*.json"],
        # "Examples": ["*.ipynb"],
    },
    include_package_data=True,
    keywords=["materials", "science", "machine", "automation", "beamline"],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    # entry_points={
    #     'console_scripts': [
    #         'meg = megnet.cli.meg:main',
    #     ]
    # }
)
