from setuptools import setup, find_packages


with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()


VERSION = '0.0.5'
DESCRIPTION = 'Stochastic and Deterministic Simulation Methods Used in Computational Biology'
HOMEPAGE = "https://github.com/LoqmanSamani/biostoch"


setup(
    name="biostoch",
    version=VERSION,
    author="Loghman Samani",
    author_email="samaniloqman91@gmail.com",
    url=HOMEPAGE,
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    include_package_data=True,
    install_requires=["numpy", "matplotlib"],
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



