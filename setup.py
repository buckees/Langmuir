"""Setup the Python Packages."""
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Geometry and Mesh", # Replace with your own username
    version="0.1.0",
    author="Wei Tian",
    author_email="buckees@gmail.com",
    description="Geometry and mesh generator for structure mesh, 1D & 2D",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com",
    packages=['Geom_Mesh',],
    # packages=setuptools.find_packages(),
    install_requires=['numpy', 'matplotlib'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
