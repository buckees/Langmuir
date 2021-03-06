"""Setup the Python Packages."""
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Langmuir Project", # Replace with your own username
    version="0.2.0",
    author="Wei Tian",
    author_email=None,
    description="Langmuir project is a comprehensive plasma model.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=None,
    # packages=['packages',],
    packages=setuptools.find_packages(),
    # the package version is specified here
    # if unsolvable errors occur, try rolling back to older versions
    install_requires=['numpy==1.19.3', 'matplotlib==3.3.0', 
                      'pandas', 'scipy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
