# bioinf-utils

`bioinf-utils` is a growing repository of lightweight and handy bioinformatics scripts designed to make routine tasks easier. Whether you're converting file formats or extracting specific genomic features, `bioinf-utils` provides the tools you need to streamline your workflows.

## Features
- **`gtf_to_bed_converter.py`**: Convert GTF files to BED format with flexible filtering and customizable output columns.

## Getting Started

### Clone the Repository
Clone the repository to your local system:
```bash
git clone https://github.com/yourusername/bioinf-utils.git
cd bioinf-utils
```

### Requirements 
bioinf-utils scripts are written in Python 3.6+ and require minimal dependencies.

Install dependencies:
```bash
pip install -r requirements.txt
```

### Examples
GTF to BED Conversion
Convert a GTF file to a BED file with default settings:
```bash
python scripts/gtf_to_bed_converter.py --input examples/example.gtf --output examples/example_output.bed --feature gene
```  

Advanced usage:
```bash
python scripts/gtf_to_bed_converter.py \
    --input examples/example.gtf \
    --output examples/example_output.bed \
    --feature exon \
    --filter_key gene_type \
    --include_patterns protein_coding \
    --bed_columns detailed \
    --name_key transcript_id

```

### Testing
Run the tests using:
```bash
python -m unittest discover -s tests -p "test_*.py"
```


## Contributing
Contributions are welcome! Please fork the repository, create a feature branch, and submit a pull request.

## License
This repository is licensed under the MIT License. See the LICENSE file for details.