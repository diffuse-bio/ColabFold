import os
import subprocess
from argparse import ArgumentParser
from google.cloud import storage


def setup_paired_msa(job_fasta, intermediate_store):
    
    base = os.path.dirname(job_fasta)
    
    # read protein ID from fasta
    with open(job_fasta, 'r') as f:
        lines = f.readlines()
        id_A = lines[0][1:-1] # ignore leading > and trailing \n
        id_B = lines[2][1:-1]
        
    
    # pull aln files from GCS
    # hardcoded GCS paths
    client = storage.Client()
    bucket = client.get_bucket('diffuse-us-central1-west1')
    blob_A = bucket.blob(f'data/msas/server_msas/intermediate_store/{id_A}.aln')
    blob_B = bucket.blob(f'data/msas/server_msas/intermediate_store/{id_B}.aln')

    assert blob_A.exists() and blob_B.exists(), f'{id_A}: {blob_A.exists()}, {id_B}: {blob_B.exists()}'
    
    # save to correct names
    with open(f'{base}/res_exp_realign.0', 'w') as f:
        f.write(blob_A.open())
    with open(f'{base}/res_exp_realign.1', 'w') as f:
        f.write(blob_B.open())
    
    # copy over res_exp_realign dbtype file
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