import setuptools

with open( "README.md", "r" ) as fh:
    long_description = fh.read()

setuptools.setup(
    name = "Cloci",
    version = "0.0b1",
    author = "xonq",
    author_email = "konkelzach@protonmail.com",
    description = "Cloci evolution-based gene cluster detection",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://gitlab.com/xonq/cloci",
    package_dir={"": "./"},
    packages = setuptools.find_packages(),
    scripts = ["cloci/cloci"], # add more
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires = '>=3.0,<4'
)
