import setuptools

with open( "README.md", "r" ) as fh:
    long_description = fh.read()

setuptools.setup(
    name = "cloci",
    version = "0.0.10",
    author = "xonq",
    author_email = "konkelzach@protonmail.com",
    description = "Function-agnostic gene cluster detection",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/xonq/cloci/",
    package_dir={"": "cloci"},
    packages = setuptools.find_packages( where="cloci" ),
    scripts = ['README.md', 'TODO.md', 'cloci/lib/output_data.py', 'cloci/lib/evo_conco.py', 'cloci/lib/generate_nulls.py', 'cloci/lib/hgp2hgx.py', 'cloci/lib/hgx2hlgs.py', 'cloci/lib/input_parsing.py', 'cloci/lib/treecalcs.py', 'README.md', 'TODO.md'],
    install_requires = ['mycotools', 'numpy', 'scipy', 'graph-tools', 'cogent3'
                        'sniffio', 'pydantic', 'PyYAML', 'SQLAlchemy', 'packaging',
                        'pydantic', 'anyio', 'attrs'],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires = '>=3.0,<4'
)
