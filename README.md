# De Bruijn Graph Assembler

## Description

This Python program is a **De Bruijn Graph Assembler** designed for reconstructing sequences from FASTQ files. The program takes sequencing reads from a FASTQ file, constructs k-mers, and builds a **De Bruijn graph** to assemble contiguous sequences (contigs). Additionally, it simplifies the graph by removing bubbles, unwanted entry and exit tips, and generates the output in **FASTA** format. It also provides an option to visualize the graph.

This program is distributed under the terms of the GNU General Public License v3.0.

## Features

- Reads sequencing data from a FASTQ file
- Constructs k-mers of customizable length
- Builds a **De Bruijn graph** from the k-mers
- Detects and simplifies graph structures:
  - Removes **bubbles** in the graph
  - Solves **entry tips** and **out tips**
- Extracts and saves **contigs** in **FASTA** format
- Visualizes the graph as an image (optional)

## Installation
### Clone the repository
Clone this repository and make sure you have Python 3 installed along with the required libraries.

```bash
git clone <repository-url>
cd <repository-directory>
```
### Install the required libraries

```
pip3 install --user networkx pytest pylint pytest-cov matplotlib
```
## Usage

Run the program with the following command:
```bash
python3 path/to/debruijn.py -i <data/input_fastq> -k <kmer_size> -o <results/output_file> [-f <results/graph_image_file>]
```
### Arguments
- `-i` : Path to the input FASTQ file (required).
- `-k` : K-mer size (default: 22).
- `-o` : Path to save the output contigs in FASTA format (default: contigs.fasta).
- `-f` : Path to save the De Bruijn graph as an image in PNG format (optional).

### Example
```bash
python3 path/to/debruijn.py -i data/sample.fastq -k 22 -o results/contigs.fasta -f results/graph.png
```

This will generate a FASTA file output.fasta containing the assembled contigs and save the graph as graph.png.

## How it works

1. K-mer Construction: The program cuts the input reads into k-mers of size k (default: 22).
2. Graph Construction: A De Bruijn graph is built where nodes represent k-1 mers, and edges represent the overlap between them.
3. Graph Simplification:
    - Bubbles are detected and removed to simplify the graph.
    - Entry and exit tips are cleaned to remove short erroneous paths.
4. Contig Extraction: Contiguous sequences (contigs) are extracted from the graph.
5. Graph Visualization: Optionally, the graph can be visualized and saved as a PNG image.

## Output

- FASTA file containing the assembled contigs.
    -Example of an output contig in FASTA format:
    ```
    >contig_0 len=1024
ATCGTACG... (contig sequence)
    ```
- Graph Image (Optional): A PNG image representing the De Bruijn graph, with large edges highlighted.

## Going further
To evaluate the quality of your assembly results (contigs), you can blast them against the reference database eva71.fna. This will help you determine how well your assembler performed. Follow the steps below:
1. Create a BLAST database from the reference file:
```bash
makeblastdb -in data/eva71.fna -dbtype nucl
```
2. Run BLAST to compare your contigs against the reference database:
```bash
blastn -query <contig_file> -db eva71.fna > output_blast.txt
```
Replace <contig_file> with the path to your assembled contig file. The results will be saved in output_blast.txt

## License
This program is licensed under the GNU General Public License v3.0. For more details, see [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.html).

## Author
- Author: Sana GUEDOUAR
- University: University of Paris-Cité
- Email: [sana.guedouar@etu.u-paris.fr](mailto:sana.guedouar@etu.u-paris.fr)