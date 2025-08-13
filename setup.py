from setuptools import setup, find_packages

setup(
    name="FMOPhore",
    version="0.1",
    author="Peter E.G.F. Ibrahim",
    author_email="pibrahim001@dundee.ac.uk, peteregfi@gmail.com",  
    description="FMOPhore for hotspot identification and classification",
    url="https://github.com/PeterEGFIbrahim/FMOPhore",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "": ["*.so"],  # Include compiled .so files
    },
    install_requires=[
        "numpy",
        "tqdm",
        "timeout-decorator",
        "argparse",
        "pandas",           
        "seaborn",
        "matplotlib",
        "rdkit"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)", 
        "Operating System :: POSIX :: Linux",
    ],
    python_requires=">=3.6",  
    entry_points={
        "console_scripts": [
            "fmophore=FMOPhore.FMOPhore_run:main",
        ],
    },
)
