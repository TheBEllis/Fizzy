#Fizzy
=======

Fizzy is a MOOSE application wrapping the nuclear inventory, source term and multi-physics code, FISPACT-II.

## Dependencies
- [FISPACT-II](https://www.ukaea.org/service/fispact/)(5.0+) binaries.
- A working [MOOSE](https://mooseframework.inl.gov/) build.
    - Being a MOOSE based application, Fizzy also has the same [minimum requirements](https://mooseframework.inl.gov/sqa/minimum_requirements.html) as MOOSE.
- PugiXML
    - PugiXML can be installed easily from most package managers
    `sudo apt install libpugixml`


## Install

To install Fizzy, users must specify the directory containg their FISPACT-II binaries. This can be done using the environment variable `FISPACT_DIR`. The exact directory will differe depending on the version of FISPACT-II being used. 

Users must also make sure that the relevent FISPACT-II libraries are present in their `LD_LIBRARY_PATH`.

#### FISPACT-II V5.1

```language=bash
export FISPACT_DIR=/path/to/FISPACT/api/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/FISPACT/api/lib/linux
```

#### FISPACT-II V5.0

```language=bash
export FISPACT_DIR=/path/to/FISPACT/ubuntu/20.10
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/FISPACT/ubuntu/20.10/lib
```

> [!IMPORTANT]
> The exact directories that need to be set may differ depending on users directory structure. By setting `FISPACT_DIR`, Fizzy will automatically look for FISPACT libraries and includes at `${FISPACT_DIR}/lib/linux` and `${FISPACT_DIR}/includes/cpp`. Users can manually set the locations of libraries and includes by settings the environment variables `FISPACT_LIB_DIR` and `FISPACT_INCLUDE_DIR`.

The final step is to run `make` in the Fizzy root directory. After a successful install, an executable called `Fizzy-opt` should be found in the application root directory. 

```language=bash
make -j 4
```


## Building docs

## 



For installation instructions and tutorials on how to use Fizzy, please visit the github pages page:


