import argparse
import os
import sys
import re
from collections import defaultdict
from src import toxinHunterHelper

th = toxinHunterHelper.ToxinHunter()

parser = argparse.ArgumentParser(description='Toxin Hunter: find potential toxin from de novo transcriptome assembly')
parser.add_argument("f", nargs="?",
                    type=argparse.FileType('r'),
                    help="Fasta file of transcriptome assembly.")
parser.add_argument("t", nargs='?',
                    type=argparse.FileType('r'),
                    help="Toxins dataset (Fasta format).")
parser.add_argument("cov", nargs='?',
                    type=float,
                    action="store",
                    help="Percentage of coverage used to select contigs in tblastn result (Value must be between 0 and 100).")
parser.add_argument("id", nargs='?',
                    action="store",
                    help="Sample name.")
args = parser.parse_args()

# Creating folders to store the results
out_dir_name = '{}_out_dir'.format(args.id)

if not os.path.exists(out_dir_name):
    out_dir_name = '{}_out_dir'.format(args.id)
    th.create_folder(out_dir_name)
    th.create_folder(os.path.join(out_dir_name,'tblastn_database'))
    th.create_folder(os.path.join(out_dir_name,'toxin_database'))

# Getting paths
path_tblastn_db = os.path.join(out_dir_name,'tblastn_database')
path_toxin_db = os.path.join(out_dir_name,'toxin_database')
path_transcripts = "{}/{}".format(path_tblastn_db, args.f.name.split("/")[-1])
path_toxins = "{}/{}".format(path_toxin_db, args.t.name.split("/")[-1])

# Copy files
os.system("cp {} {}" .format(args.f.name, path_tblastn_db))
os.system("cp {} {}" .format(args.t.name, path_toxin_db))

# Formatting blast database
th.format_blast_nts_db(path_transcripts)
th.format_blast_pts_db(path_toxins)

# tBlastn
th.run_blast('tblastn', path_transcripts, path_toxins, out_dir_name, args.id)
print("\n\n*** Step 1: tblastn - ok ***\n\n")

# Getting toxin list
th.get_toxin_list(out_dir_name, args.id)

# Filtering all blast hits and hsps from tblastn
toxin_list = "{}/list_toxins".format(out_dir_name)
tblastn_result = "{}/{}.tblastn.out".format(out_dir_name, args.id)
result_filter = "{}/{}_filtered_tblastn.tab".format(out_dir_name, args.id)
th.filter_all_hits_hsps_from_tblastn(toxin_list, tblastn_result, args.cov, result_filter)

# Getting selected transcripts
th.get_selected_transcripts(out_dir_name, args.id)
transcripts_list = "{}/selected_transcripts_from_tblastn".format(out_dir_name)
output_selected_fasta = "{}/{}_selected_contigs_from_tblastn.fasta".format(out_dir_name, args.id)
th.get_fasta_selected_transcripts(transcripts_list, path_transcripts, output_selected_fasta)
print("\n\n*** Step 2: Select transcripts of potential toxins from whole transcriptome - ok ***\n\n")

# Getting ORFs from selected transcripts
selected_fasta = output_selected_fasta.split("/")[-1]
os.chdir(out_dir_name)
th.run_transdecoder(selected_fasta)
os.chdir('..')

print("\n\n*** Step 3: Get ORFs from selected transcripts - ok ***\n\n")

# Blastp: selectin and annotating the proteins
blastp_query_path = '{}/{}_selected_contigs_from_tblastn.fasta.transdecoder_dir/longest_orfs.pep'.format(out_dir_name, args.id)
th.run_blast('blastp', path_toxins, blastp_query_path, out_dir_name, args.id)
print("\n\n*** Step 4: Blastp - ok ***\n\n")


