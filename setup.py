from setuptools import setup, find_packages

setup(
    name="FMOPhore",
    version="0.1",
    author="Peter E.G.F. Ibrahim",
    author_email="2448959@dundee.ac.uk peteregfi@gmail.com",
    description="FMOPhore for hotspot identification and classification",
    url="https://github.com/PeterEGFIbrahim/FMOPhore.git", 
    packages=find_packages(),
    install_requires=[
        "numpy",
        "tqdm",
        "timeout-decorator",
        "argparse",
    ], 
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0 license",
        "Operating System :: Linux",
    ],
    python_requires=">=3.6", 
    entry_points={
        "console_scripts": [
            "fmophore=FMOPhore.FMOPhore_run:main",
        ],
    },
)
