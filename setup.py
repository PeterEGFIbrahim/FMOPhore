from setuptools import setup, find_packages

# Optional: Read from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="FMOPhore",
    version="0.1",
    author="Peter E.G.F. Ibrahim",
    author_email="2448959@dundee.ac.uk, peteregfi@gmail.com",
    description="FMOPhore for hotspot identification and classification",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PeterEGFIbrahim/FMOPhore",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "tqdm",
        "timeout-decorator",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: Linux",
    ],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "fmophore=FMOPhore.FMOPhore_run:main",
        ],
    },
)
