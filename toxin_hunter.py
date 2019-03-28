import argparse
import os
import sys
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
print("\n*** Step 1: tblastn - ok ***\n")

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
print("\n*** Step 2: Select transcripts of potential toxins from whole transcriptome - ok ***\n")

# Getting ORFs from selected transcripts
selected_fasta = output_selected_fasta.split("/")[-1]
os.chdir(out_dir_name)
th.run_transdecoder(selected_fasta)
os.chdir('..')
print("\n*** Step 3: Get ORFs from selected transcripts - ok ***\n")

# Blastp: selecting and annotating the proteins
blastp_query_path = '{}/{}_selected_contigs_from_tblastn.fasta.transdecoder_dir/longest_orfs.pep'.format(out_dir_name, args.id)
th.run_blast('blastp', path_toxins, blastp_query_path, out_dir_name, args.id)
print("\n*** Step 4: Blastp - ok ***\n")

blastp_output = "{}/{}.blastp.out".format(out_dir_name, args.id)
output_path = "{}/{}".format(out_dir_name, args.id)
th.select_annotate_proteins(blastp_output, blastp_query_path, output_path)
print("\n*** Step 5: Filtering potential toxins based on coverage and identity - ok ***\n")

# Select nucleotides ORFs
th.get_ids_nts_orfs(out_dir_name, args.id)
th.format_nts_orfs(out_dir_name, args.id)
ids_nts_path = "{}/ids_get_nts_orfs".format(out_dir_name)
orfs_nts_path = '{}/{}_selected_contigs_from_tblastn.fasta.transdecoder_dir/longest_orfs.cds'.format(out_dir_name, args.id)
th.get_nts_orfs(ids_nts_path, orfs_nts_path, out_dir_name, args.id)
print("\n*** Step 6: Get nucleotide ORFs - ok ***\n")

# Remove repeated sequences
toxin_nts_path = "{}/{}_selected_orfs_annotated_nts.fasta".format(out_dir_name, args.id)
toxin_pts_path = "{}/{}_selected_orfs_annotated_pts.fasta".format(out_dir_name, args.id)
th.remove_redundancy(toxin_nts_path, toxin_pts_path)
print("\n*** Step 7: Remove redundancy - ok ***\n")

# Selecting whole contigs from selected ORFs
th.get_ids_whole_contigs(out_dir_name, args.id)
th.get_whole_contigs_from_selected_orfs(out_dir_name, args.id)
print("\n*** Step 8: Get whole contigs from selected potential toxins - ok ***\n")

# Counting toxin families in the final results
fasta_file_path = "{}/{}_selected_orfs_annotated_nts_uniq.fasta".format(out_dir_name, args.id)
th.count_toxin_families(fasta_file_path, out_dir_name, args.id)
print("\n*** Step 9: Counting toxin families - ok ***\n")

# Removing intermediary files
os.chdir(out_dir_name)
os.system("rm -r tblastn_database toxin_database list_toxins selected_transcripts_from_tblastn ids_get_nts_orfs \
           *_selected_orfs_annotated_pts.fasta *_selected_orfs_annotated_nts.fasta ids_get_whole_contigs")
print("\n*** Step 10: Removing intermediary files - ok ***\n")
print("*** FINISH SUCCESSFULLY!***\n")
