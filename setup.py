import setuptools

setuptools.setup(
    name="pypart",
    version="1.0.0",
    author="Ze Zhao",
    author_email="zezhao2@illinois.edu",
    description="A Python package for mesh partition",
    packages=setuptools.find_packages(),
    install_requires=[
        "numpy",
        "meshio",
        "h5py",
        "pymetis>=2023.1",
    ],
    entry_points={
        "console_scripts": [
            "pypart=pypart.pypart:main",
        ],
    },
)

