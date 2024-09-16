import os
import subprocess
from argparse import ArgumentParser
from google.cloud import storage
from add_tax_to_msa import parse_fasta


def setup_single_msa(job_fasta):
    
    base = os.path.dirname(job_fasta)
    ids = parse_fasta(job_fasta)[1]

    # pull aln files from GCS; hardcoded GCS paths
    client = storage.Client()
    bucket = client.get_bucket('diffuse-us-central1-west1')
    blobs = []
    for i in ids:
        blobs += [bucket.blob(f'data/msas/server_msas/intermediate_store/{i}.aln')]
       
    # save file from gcs to pick up where we left off 
    # also dynamically build the index as the files are copied
    index_str = ''
    total_size = 0
    # check that all files exist
    exists = [blob.exists() for blob in blobs]
    print(f'Out of {len(ids)} IDs, aln files were found for {sum(exists)} of them.')
    
    # only copy over intermediate files to local if running the second step
    if all(exists):
        subprocess.run(['touch', f'{base}/ALN_FOUND'])
        for j, blob in enumerate(blobs):
            with open(f'{base}/res_exp_realign.{j}', 'w') as f:
                f.write(blob.open('r').read())
            # build index string
            curr_size = os.stat(f'{base}/res_exp_realign.{j}').st_size
            index_str += f'{j}\t{total_size}\t{curr_size}\n'
            total_size += curr_size   

        with open(f'{base}/res_exp_realign.index', 'w') as f:
            f.write(index_str)

        subprocess.run(['cp', f'{os.getcwd()}/res_exp_realign.dbtype', f'{base}/'])
    else:
        print('Calculating intermediate files...')

def main():
    parser = ArgumentParser()
    parser.add_argument(
        "job_fasta",
        type=str,
        help="Path to fasta file that contains the protein",
    )
    args = parser.parse_args()
    setup_single_msa(args.job_fasta)


if __name__ == '__main__':
    main()
