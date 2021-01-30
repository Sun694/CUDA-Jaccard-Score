import rdkit
from rdkit import Chem
import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input file path", required=True)
parser.add_argument("-o", default="sanitized_smiles.csv", help="Output file path", required=False)
args = vars(parser.parse_args())



def canon(x):
    x = Chem.MolFromSmiles(x)
    if x is None:
        return x
    else:
        return Chem.MolToSmiles(x)
        

if __name__ == "__main__":
    file = args["i"]
    data = pd.read_csv(file, header=None)[:30]
    data.head()
    data = data.dropna()
    # some invalid SMILES have spaces in them
    data = data[0].map(lambda x: np.nan if " " in x else x)
    # Drop indices where the SMILES string has a space in it
    data = data.dropna()
    data = data.reset_index(drop=True)
    # Drop indices where the target/molecule pair are duplicates
    data = data.drop_duplicates()
    data = data.reset_index(drop=True)
    print(data)
    # Drop indices where the SMILES string is invalid
    data = data.apply(lambda x: canon(x))
    data = data.dropna()
    data = data.reset_index(drop=True)
    data.to_csv(args["o"], header=None, index=False)