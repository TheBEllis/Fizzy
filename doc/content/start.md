# Installation

## Prerequisites 

The 2 things you will need to use Fizzy are:

1. A +FISPACT-II+ installation
2. A working +MOOSE+ installation

### Optional 

To use Fizzy's distributed sampling capabilities, users will also need a valid Cardinal installation.

## Installing

Installing Fizzy is very similar to installing any other MOOSE based application. The only extra requirements are some environement variables.

First at build time, the environment variable `FISPACT_DIR` must be set to the directory containing a valid set of FISPACT libraries and includes. For example:

```language=bash
export FISPACT_DIR=../../FISPACT/ubuntu/20.10
```

At run time, the relevent FISPACT-II libraries must be present in your `LD_LIBRARY_PATH`.

```language=bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../FISPACT/ubuntu/20.10/lib
```
