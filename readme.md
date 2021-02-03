# FastJacc
 A fast molecular database search function that offers over a 100x times speedup over RDKit for bulk tanimoto.

Using an optimized custom CUDA kernel, on a P100, you can compute tanimoto between 800,000 molecules in less than 4ms.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites

Only dependencies are numpy, rdkit, gcc, and CUDA >= 10.1.

`pip install {numpy, pandas}`

`conda install -c rdkit rdkit`

For gcc / CUDA go to the respective websites online.

### Usage

Sanatize.py: 
```
usage: sanitize.py [-h] -i INPUT_DIRECTORY // input CSV for process. assumes it is comma delimited, one molecule per line.
                     -o OUTPUT_DIRECTORY // outputs in the same format as above

Sanitize a given input CSV, destroying all entries that RDKit cannot process.
```

convert_to_binary.py:
```
usage: convert_to_binary.py [-h] -i INPUT_DIRECTORY
                     -o OUTPUT_DIRECTORY
                     -n CHUNK_SIZE
Convert the given sanitized input into a binary format for fast processing.
```

fastSearch_CUDA.cu:

```
usage: fastSearch_CUDA DATABASE_DIRECTORY QUERY_DIRECTORY OUTPUT_DIRECTORY BLOCK_SIZE TOP_K

Query the database with the given query vectors, gather the top k similarities, then output them in human-readable format in output directory.

```

### Example usage:

```
python sanitize.py -i database_smiles_sequences.csv -o sanitized_database.csv
python convert_to_binary.py -i sanitized_database.txt -o database.bin 
python sanitize.py -i query_smiles_sequences.csv -o sanitized_queries.csv
python convert_to_binary.py -i sanitized_queries.txt -o queries.bin 

nvcc fastSearch_CUDA.cu -o fastSearch

fastSearch database.bin queries.bin search_results.txt 1024 30
```

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
