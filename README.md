# cube_processor

A utility script for working with Gaussian cube files.

## Requirements

- Python 3
- Numpy

## Usage

Up-to-date usage instructions can always be obtained by executing the script with the `--help` option.
```
usage: cube_processor.py [-h] [--print-info] [--calc-iso-value [PERCENTAGE]] cube_file

CLI tool for processing Gaussian cube files

positional arguments:
  cube_file             The cube file to process

options:
  -h, --help            show this help message and exit
  --print-info          Print some general info about the read cube file
  --calc-iso-value [PERCENTAGE]
                        Determine a threshold value for a n isosurface such that the enclosed volume will contain the given percentage of the total property
```

