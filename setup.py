from setuptools import setup, find_packages

VERSION = "3.0"
DESCRIPTION = "Cellori"
LONG_DESCRIPTION = "A fast and robust algorithm for clustered nuclei segmentation."

setup(
    name="cellori",
    version=VERSION,
    author="William Niu",
    author_email="<wniu721@gmail.com>",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=['matplotlib', 'numba', 'numpy', 'opencv-python', 'pyside6', 'scikit-image', 'scipy', 'simpleitk',
                      'stitchwell', 'tifffile'],
    keywords=["nuclei", "segmentation"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)
