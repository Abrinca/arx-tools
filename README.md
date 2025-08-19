# Arx Tools

A set of scripts that helps to import genome data into the Arx folder structure.

## Quick Start

```bash
# Install
pip install arx-tools

# Initialize folder structure
export FOLDER_STRUCTURE=/path/to/folder_structure
init_folder_structure

# Import genome data
import_genome --import_dir=/path/to/genome/data --organism STRAIN --genome STRAIN.1
```

## Documentation

📖 **Documentation:** [docs.abrinca.com/tools](https://docs.abrinca.com/tools/)

## Help

All scripts have built-in help functions:

```bash
import_genome --help
```

## Requirements

- Python 3.9 or higher
