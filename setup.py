from setuptools import setup, find_packages




VERSION = '0.0.1'
DESCRIPTION = 'Stochastic and Deterministic Simulation Methods Used in Computational Biology'
LONG_DESCRIPTION = README.md



HOMEPAGE = "https://github.com/LoqmanSamani/biostoch"

# Setting up
setup(
    name="BioStoch",
    version=VERSION,
    author="Loghman Samani",
    author_email="samaniloqman91@gmail.com",
    url=HOMEPAGE,
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=["pandas", "numpy", "seaborn", "matplotlib"],
    keywords=['Stochastic Simulation', 'Deterministic Simulation', 'Computational Biology'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)



