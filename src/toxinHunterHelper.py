import os
from Bio import SeqIO

class ToxinHunter:
  # def __init__(self, folder_name):
  #     self.folder_name = folder_name

  def create_folder(self, folder_name):
    os.makedirs(folder_name)

  def format_blast_nts_db(self, database):
    os.system("makeblastdb -in {} -parse_seqids -dbtype nucl".format(database))

  def format_blast_pts_db(self, database):
    os.system("makeblastdb -in {} -parse_seqids -dbtype prot".format(database))

  def run_blast(self, blast_program, database, query, out_dir_name, sample_name):
    os.system("{} -db {} -query {} -evalue 1e-6 -out {}/{}.tblastn.out -num_threads 20 \
          -outfmt '6 qseqid qlen qframe sseqid slen sframe pident length mismatch gapopen qstart qend sstart send qcovs qcovhsp evalue bitscore'\
          ".format(blast_program, database, query, out_dir_name, sample_name))

  def get_toxin_list(self, out_dir_name, sample_name):
    os.system("cat {}/{}.tblastn.out | cut -f 1 | sort | uniq > {}/list_toxins"\
        .format(out_dir_name, sample_name, out_dir_name))

  def filter_all_hits_hsps_from_tblastn(self, path_toxin_list, path_tblastn_result, coverage, path_result_filter):
    toxin_list = open(path_toxin_list, 'r')
    tblastn_result = path_tblastn_result
    cov_filter = coverage
    result_tblastn_filter = open(path_result_filter, "w")

    dict_toxin = {}
    for t in toxin_list:
        toxin = t.strip()
        blast_reports = open(tblastn_result, 'r')
        dict_toxin[toxin] = []
        for r in blast_reports:
            report = r.strip()
            toxin_in_tab = report.split('\t')[0]

            if toxin.strip() == toxin_in_tab:
                dict_toxin[toxin].append(report)
        blast_reports.close()

    final_dict_trans_by_toxin = {k:v for k,v in dict_toxin.items() if v}

    # Selecting list of hits
    my_final_dict = {}

    for k,v in final_dict_trans_by_toxin.items():
        hits_list = {}
        for i in v:
            hit_fields = i.split('\t')
            qlen = hit_fields[1]
            qstart = hit_fields[10]
            qend = hit_fields[11]
            coord = ' '.join([qstart, qend])

            if hit_fields[3] in hits_list.keys():
                hits_list[hit_fields[3]].append(coord)
            else:
                hits_list[hit_fields[3]] = [coord]

        for keys, values in hits_list.items():
            name = keys + "|" + k + "|" + qlen
            my_final_dict[name] = values

    for i, j in my_final_dict.items():
        len_query = i.split("|")[-1]
        query = i.split("|")[1]
        transcript = i.split("|")[0]

        #Merging intervals and calculating coverage
        if len(j) > 1:
            coord_list = []
            for pair in j:
                start, end = pair.split(" ")
                pair_in_list = [int(start),int(end)]
                coord_list.append(pair_in_list)

            sorted_intervals = sorted(coord_list, key=lambda x: x[0])
            merged =[]

            for x_list in sorted_intervals:
                if not merged:
                    merged.append(x_list)
                else:
                    b = merged.pop()
                    if int(b[1]) >= int(x_list[0]) or (int(b[1]) + 1) ==  int(x_list[0]):
                        new_x_list = list([b[0], x_list[1]])
                        merged.append(new_x_list)
                    else:
                        merged.append(b)
                        merged.append(x_list)
            cov = 0.0
            if len(merged) > 1:
                len_alig_total = 0

                for each_pair in merged:
                    len_alig = (int(each_pair[1]) - int(each_pair[0])) + 1
                    len_alig_total += len_alig
                cov = (int(len_alig_total) / int(len_query)) *100

            else:
              len_alig = (merged[0][1] - merged[0][0]) + 1
              cov = (int(len_alig) / int(len_query)) *100

            if cov >= float(cov_filter):
                result_tblastn_filter.write("{}\t{}\t{}\t{}\t{}\n" .format(query,transcript,merged,len_query,cov))
        else:
            start, end = j[0].split(" ")
            len_alig_total = (int(end) - int(start)) + 1
            cov = (int(len_alig_total) / int(len_query)) *100

            if cov >= float(cov_filter):
                result_tblastn_filter.write("{}\t{}\t{}\t{}\t{}\n" .format(query,transcript,j,len_query,cov))

    toxin_list.close()
    result_tblastn_filter.close()

  def get_selected_transcripts(self, out_dir_name, sample_name):
    os.system("cat {}/{}_filtered_tblastn.tab | cut -f 2 | sort | uniq > {}/selected_transcripts_from_tblastn"\
        .format(out_dir_name, sample_name, out_dir_name))

  def get_fasta_selected_transcripts(self, list_selected_transcripts, fasta_all_transcripts, path_result):
    transcripts_list = open(list_selected_transcripts ,'r')
    outfile = open(path_result, "w")

    mydict = {}

    for i in transcripts_list:
      list_item = i.strip()
      input_fasta = open(fasta_all_transcripts, "r")
      fasta_sequences = SeqIO.parse(input_fasta,'fasta')

      for fasta in fasta_sequences:
          name, sequence = fasta.id, str(fasta.seq)
          seq_id = name.split(" ")[0]

          if list_item == seq_id:
              mydict[seq_id]= sequence
      input_fasta.close()

    for k,v in mydict.items():
      outfile.write(">{}\n{}\n".format(k,v))

    transcripts_list.close()
    outfile.close()

  def run_transdecoder(self, fasta_file):
    print("TransDecoder.LongOrfs -t {} -m 50".format(fasta_file))
    os.system("TransDecoder.LongOrfs -t {} -m 50".format(fasta_file))
    os.system("rm -r *.cmds *.transdecoder_dir.__checkpoints_longorfs")





