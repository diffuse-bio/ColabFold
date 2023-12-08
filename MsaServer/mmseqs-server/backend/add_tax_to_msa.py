import pandas as pd
import os
from argparse import ArgumentParser
import glob

def parse_fasta(fasta_string: str):
    """Parses FASTA string and returns list of strings with amino-acid sequences.

    Arguments:
      fasta_string: The string contents of a FASTA file.

    Returns:
      A tuple of two lists:
      * A list of sequences.
      * A list of sequence descriptions taken from the comment lines. In the
        same order as the sequences.
    """
    if os.path.exists(fasta_string):
        fasta_string = open(fasta_string,'r').read()
    
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append("")
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line

    return sequences, descriptions


def parse_fasta_as_df(fasta_string):
    '''
    Instead of returning the fasta as two lists, convert it to a dataframe
    '''
    msa = parse_fasta(fasta_string)
    sep_char = '\t'
    if '|' in msa[1][-1]:
        sep_char = '|'
    rows = [[msa[0][i]]+msa[1][i].split(sep_char) for i in range(len(msa[0]))]
    msa_df = pd.DataFrame(rows)
    return msa_df


def to_fasta(seq_list, seq_id_list, output_fn):
    
    # add '>' to headers if not already prepended
    seq_id_list = [f'>{seq_id}' if '>'!=seq_id[0] else seq_id for seq_id in seq_id_list]
    newlines = [val for pair in zip(seq_id_list, seq_list) for val in pair]
    with open(output_fn, 'w') as f:
        f.write('\n'.join(newlines))
    return


def add_tax_to_msa(convertalis_path: str, msa_path: str):
    '''
    Appends taxonomy ID, name, and lineage from `mmseqs convertalis --format-output target,evalue,taxid,taxname,taxlineage`
    to the msa file(s), by matching fasta IDs. 
    If the msa_path is a directory, appends taxonomy info to all .a3m files found.
    Also converts headers to diffpalm format (expects species to be in header.split('|')[1])
    '''
    convertalis = pd.read_csv(convertalis_path, sep='\t', header=None)
    convertalis.columns = ['target', 'evalue', 'taxid', 'taxname', 'taxlineage']
    convertalis['taxname'] = convertalis['taxname'].str.replace(' ','_')
    convertalis['taxlineage'] = convertalis['taxlineage'].str.replace(' ','_')
    convertalis = convertalis.drop_duplicates()
    
    if os.path.isdir(msa_path):
        fns = glob.glob(f'{msa_path}/*.a3m') # batch add tax for all a3ms found
    else:
        fns = [msa_path]

    for fn in fns:
        with open(fn,'r') as f:
            msa_df = parse_fasta_as_df(f.read())
            msa_df = msa_df.rename(columns={0:'seq', 1:'id'})
            
            # merge convertalis and a3m file
            m = msa_df.merge(convertalis, how='left', left_on='id', right_on='target').astype(str)
            m = m.drop(columns='target')
            out_a3m0 = list(m['seq'])
            out_columns = ['id','taxname','taxid','evalue'] + list(msa_df.columns)[2:] + ['taxlineage']
            out_a3m1 = list(m[out_columns].agg('|'.join, axis=1))
            to_fasta(out_a3m0, out_a3m1, f'{fn}.tax')
    return m


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "convertalis_path",
        type=str,
        help="Path to convertalis_tax file",
    )
    parser.add_argument(
        "msa_path",
        type=str,
        help="Path to MSA to append tax info to",
    )
    args = parser.parse_args()
    add_tax_to_msa(args.convertalis_path, args.msa_path)


if __name__ == '__main__':
    main()
