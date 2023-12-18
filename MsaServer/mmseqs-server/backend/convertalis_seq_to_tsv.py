import os
import subprocess
from argparse import ArgumentParser
from google.cloud import storage
import pandas as pd

def convertalis_seq_to_tsv(convertalis_seq, a3m_0, a3m_1):
    with open(convertalis_seq, 'r') as convertalis:
        lines = convertalis.readlines()
        seq_dict = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in lines}
        
        # read in the two a3m files
        a3m_0 = open(a3m_0, 'r').readlines()[0::2]
        a3m_0 = [line.strip().split('>')[1] for line in a3m_0]

        a3m_1 = open(a3m_1, 'r').readlines()[0::2]
        a3m_1 = [line.strip().split('>')[1] for line in a3m_1]
        
        seq_df = pd.DataFrame(seq_dict.items(), columns=['id', 'seq'])
        new_df = pd.DataFrame([a3m_0, a3m_1]).T
        new_df.columns = ['id_1', 'id_2']
        
        # merge dataframes to map id to seq
        merge_df = new_df.merge(seq_df, left_on='id_1', right_on='id').merge(seq_df, left_on='id_2', right_on='id')
        merge_df = merge_df[['id_1', 'seq_x', 'id_2', 'seq_y']]
        merge_df.columns = ['id_1', 'seq_1', 'id_2', 'seq_2']

        merge_df.to_csv(f'{convertalis_seq}.tsv', sep='\t', index=False)


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "convertalis_seq",
        type=str,
        help="Path to convertalis_seq that contains the two proteins in the PPI",
    )
    parser.add_argument(
        "a3m_0",
        type=str,
        help="Path to 0.a3m",
    )
    parser.add_argument(
        "a3m_1",
        type=str,
        help="Path to 1.a3m",
    )
    args = parser.parse_args()
    convertalis_seq_to_tsv(args.convertalis_seq, args.a3m_0, args.a3m_1)


if __name__ == '__main__':
    main()
