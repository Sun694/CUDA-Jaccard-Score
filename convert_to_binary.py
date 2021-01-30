import rdkit
from rdkit import Chem
import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input file path. Heavily recommend sanitizing via sanitize.py first", required=True)
parser.add_argument("-o", default="mols.bin", help="Output file path", required=False)
parser.add_argument("-n", default=32, help="Number of bits in a chunk. Must be a factor of 2 less than 64. Modifying this requires modifying corresponding C / CUDA code manually. Recommended do not touch.", required=False)

args = vars(parser.parse_args())

def fingerprintify(x):
    return Chem.RDKFingerprint(Chem.MolFromSmiles(x))
    
def bitstringify(x):
    return x.ToBitString()
    
def int2bytes(x):
    return (x).to_bytes(length=4, byteorder="little")
    
if __name__ == "__main__":
    file = args["i"];
    data = pd.read_csv(file, header=None)[0]
    print(data)
    data = data.apply(lambda x: fingerprintify(x))
    data = data.apply(lambda x: bitstringify(x))
    with open(args["o"], "wb") as f:
        f.write(int2bytes(len(data)))
        for i in range(len(data)):
            current_bit_vect = data.iloc[i]
            num_on = current_bit_vect.count("1")
            # bits per chunk
            n = args["n"]
            # breaks bitstring into chunks that can be converted to ints
            chunked_vec = [current_bit_vect[idx : idx + n] for idx in range(0, len(current_bit_vect), n)]
            chunked_vec = [int(x, 2) for x in chunked_vec]

            f.write(int2bytes(num_on))
            [f.write(int2bytes(x)) for x in chunked_vec]
