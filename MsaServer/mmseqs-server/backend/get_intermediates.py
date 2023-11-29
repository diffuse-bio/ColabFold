import os
import subprocess
from argparse import ArgumentParser


def setup_paired_msa(job_fasta, intermediate_store):
    
    base = os.path.dirname(job_fasta)
    
    # read protein ID from fasta
    with open(job_fasta, 'r') as f:
        lines = f.readlines()
        id_A = lines[0][1:-1] # ignore leading > and trailing \n
        id_B = lines[2][1:-1]
    

    # copy over res_exp_realign files -- symbolic link?
    subprocess.run(['ln', '-s', f'{intermediate_store}/{id_A}.aln', f'{base}/res_exp_realign.0'])
    subprocess.run(['ln', '-s', f'{intermediate_store}/{id_B}.aln', f'{base}/res_exp_realign.1'])
    subprocess.run(['ln', '-s', f'{intermediate_store}/res_exp_realign.dbtype', f'{base}/'])



    # remake the index
    index_a = os.stat(f'{base}/res_exp_realign.0').st_size
    index_b = os.stat(f'{base}/res_exp_realign.1').st_size
    with open(f'{base}/res_exp_realign.index', 'w') as f:
        f.write(f'0\t0\t{index_a}\n')
        f.write(f'1\t{index_a}\t{index_b}\n')


def main():
    parser = ArgumentParser()
    parser.add_argument(
        "job_fasta",
        type=str,
        help="Path to fasta file that contains the two proteins in the PPI",
    )
    parser.add_argument(
        "intermediate_store",
        type=str,
        help="Path to directory that contains all intermediate files, each named ID.aln",
    )
    args = parser.parse_args()
    setup_paired_msa(args.job_fasta, args.intermediate_store)



if __name__ == '__main__':
    main()