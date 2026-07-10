# Fizzy

Fizzy is a MOOSE application wrapping the nuclear inventory, source term and multi-physics code, FISPACT-II.

## Dependencies
- [FISPACT-II](https://www.ukaea.org/service/fispact/)(5.0+) binaries.
- A working [MOOSE](https://mooseframework.inl.gov/) build.
    - Being a MOOSE based application, Fizzy also has the same [minimum requirements](https://mooseframework.inl.gov/sqa/minimum_requirements.html) as MOOSE.
- PugiXML
    - PugiXML can be installed easily from most package managers
    - `sudo apt install libpugixml`

## Obtaining Fizzy

To obtain Fizzy, simply clone this repository.

```language=bash
git clone https://github.com/TheBEllis/Fizzy.git
cd Fizzy
```

## Install

> [!WARNING]
> Fizzy has only been tested on Ubuntu 24.04 (so far). Other Linux distros, MacOS and Windows have not been tested. If users trying to run Fizzy on these platforms
experience issues with installation and/or running, please report the issue on Fizzy's [github issues](https://github.com/TheBEllis/Fizzy/issues) page.

To install Fizzy, users must specify the directory containg their FISPACT-II binaries. This can be done using the environment variable `FISPACT_DIR`. The exact directory will differe depending on the version of FISPACT-II being used. 

Users must also make sure that the relevent FISPACT-II libraries are present in their `LD_LIBRARY_PATH`.

#### For FISPACT-II V5.1

```language=bash
export FISPACT_DIR=/path/to/FISPACT/api/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/FISPACT/api/lib/linux
```

#### For FISPACT-II V5.0

```language=bash
export FISPACT_DIR=/path/to/FISPACT/ubuntu/20.10
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/FISPACT/ubuntu/20.10/lib
```

> [!IMPORTANT]
> The exact directories that need to be set may differ depending on users directory structure. By setting `FISPACT_DIR`, Fizzy will automatically look for FISPACT libraries and includes at `${FISPACT_DIR}/lib/linux` and `${FISPACT_DIR}/includes/cpp`. Users can manually set the locations of libraries and includes by settings the environment variables `FISPACT_LIB_DIR` and `FISPACT_INCLUDE_DIR`.

The final step is to run `make` in the Fizzy root directory. After a successful install, an executable called `Fizzy-opt` should be found in the application root directory. 

```language=bash
cd /path/to/Fizzy
make -j 4
```

## Building docs
Due to Fizzy currently being private, the only way to access documentation is by building it locally. Users can do this by running.

```
cd /path/to/Fizzy/doc
./moosedocs.py build --destination ./build_folder
```

This will populate the specified `build_folder` with html files. To access the documentation after building, `cd` into the build folder and open `index.html` from your chosen web browser.

## References



