import os
import subprocess
from argparse import ArgumentParser
from google.cloud import storage
import pandas as pd
import glob

def pair_seq_to_tsv(convertalis_seq, a3m_0, a3m_1):
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


def single_seq_to_tsv(convertalis_seq, a3m_dir):
    print('single_seq_to_tsv')
    fns = glob.glob(f'{a3m_dir}/*.a3m')
    with open(convertalis_seq, 'r') as convertalis:
        lines = convertalis.readlines()
        seq_dict = {line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in lines}
        seq_df = pd.DataFrame(seq_dict.items(), columns=['id', 'seq'])

        for fn in fns:
            print(fn)
            protein_id = fn.split('/')[-1].split('.')[0]
            a3m = open(fn, 'r').readlines()[0::2]
            a3m = [line.strip().split('>')[1].split('\t')[0] for line in a3m]            
            a3m_df = pd.DataFrame({'id': a3m})
            
            # merge dataframes to map id to seq
            merge_df = a3m_df.merge(seq_df, left_on='id', right_on='id')
            merge_df = merge_df[['id', 'seq']]
            merge_df.to_csv(f'{a3m_dir}/{protein_id}.tsv', sep='\t', index=False)



def main():
    parser = ArgumentParser()
    parser.add_argument(
        "convertalis_seq",
        type=str,
        help="Path to convertalis_seq that contains the two proteins in the PPI",
    )
    parser.add_argument(
        "--pair",
        action='store_true',
        help="Paired tsv format",
    )
    parser.add_argument(
        "--a3m_dir",
        type=str,
        help="Path to directory containing a3ms (named ID.a3m)",
    )
    parser.add_argument(
        "--a3m_0",
        type=str,
        help="Path to 0.a3m",
    )
    parser.add_argument(
        "--a3m_1",
        type=str,
        help="Path to 1.a3m",
    )
    args = parser.parse_args()
    if args.pair:
        pair_seq_to_tsv(args.convertalis_seq, args.a3m_0, args.a3m_1)
    else:
        single_seq_to_tsv(args.convertalis_seq, args.a3m_dir)


if __name__ == '__main__':
    main()
