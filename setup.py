import setuptools

with open( "README.md", "r" ) as fh:
    long_description = fh.read()

setuptools.setup(
    name = "bigmcl",
    version = "0.2b2",
    author = "xonq",
    author_email = "konkelzach@protonmail.com",
    description = "Large scale Markov clustering (MCL) via subgraph extraction",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://gitlab.com/xonq/bigmcl",
    package_dir={"": "./"},
    packages = setuptools.find_packages(),
    scripts = ["bigmcl/bigmcl.py"],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires = '>=3.0,<4'
)
