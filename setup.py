import numpy
from skbuild import setup

numpy_include_dir = numpy.get_include()

setup(
    name="crpropa",
    version="3.2.0",
    description="CRPropa: A Monte Carlo simulation framework for high-energy cosmic ray propagation",
    author="CRPropa Collaboration",
    url="https://github.com/CRPropa/CRPropa3",
    install_requires=["numpy"],
    package_dir={"": "python"},
    packages=["crpropa"],
    include_package_data=True,
    python_requires=">=3.7",
    # Pass the NumPy include directory (using the name expected by CMake)
    cmake_args=[f"-DPython_NumPy_INCLUDE_DIRS={numpy_include_dir}", "-DPython_FIND_NUMPY=ON", "-DDOWNLOAD_DATA=OFF",],
)