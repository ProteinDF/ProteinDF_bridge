# ProteinDF_bridge -- bridge scripts the ProteinDF package and other data/package

## What you need to use ProteinDF_bridge

* Python 2.x or later (recommended: above 3.6)
* Python modules
  * setuptools
  * numpy
  * argparse
  * PyYaml
  * msgpack-python, u-msgpack-python or msgpack-pure
  * configparser
  * ordereddict (for Python 2.6)

## How to install ProteinDF_bridge

### get source files

Clone a copy of the main ProteinDF_bridge git repo by running:

```bash
git clone git://github.com/ProteinDF/ProteinDF_bridge.git
```


### install ProteinDF_bridge module using the standard `venv` module (recommended)

- required above Python 3.3

#### confirm install directory

Here, we will assume that the virtual environment of Python will be installed in the directory specified by the environment variable `PDF_HOME`.
Set the environment variables according to the shell you are using.

The following shows how to set them in bash:

```
$ export PDF_HOME=${HOME}/local/ProteinDF
```

#### prepare standard virtual environment of Python

The following shows how to set them in bash:

```
$ python -m venv ${PDF_HOME}
$ source ${PDF_HOME}/bin/activate
```

#### install module

Enter the source directory and run build script:

```bash
$ cd ProteinDF_bridge
$ pip install --use-feature=in-tree-build .
```

# Documents

sorry, in preparation.


# License

ProteinDF_bridge is licensed under the GNU GPL v3.
The source code can be found on the Github.


# Bugs

If you find any bugs, please let me know.
And if you have a suggestion for improvement, please let me know.
