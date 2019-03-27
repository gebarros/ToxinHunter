import os
import re
from Bio import SeqIO
from collections import defaultdict

class ToxinHunter:

  def create_folder(self, folder_name):
    os.makedirs(folder_name)

  def format_blast_nts_db(self, database):
    os.system("makeblastdb -in {} -parse_seqids -dbtype nucl".format(database))

  def format_blast_pts_db(self, database):
    os.system("makeblastdb -in {} -parse_seqids -dbtype prot".format(database))

  def run_blast(self, blast_program, database, query, out_dir_name, sample_name):
    os.system("{} -db {} -query {} -evalue 1e-6 -out {}/{}.{}.out -num_threads 20 \
          -outfmt '6 qseqid qlen qframe sseqid slen sframe pident length mismatch gapopen qstart qend sstart send qcovs qcovhsp evalue bitscore'\
          ".format(blast_program, database, query, out_dir_name, sample_name,blast_program))

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

  def select_annotate_proteins(self, blastp_output, orfs_file, output_path):
    blastp_tab = open(blastp_output, 'r')
    orfs_fasta = orfs_file
    cov_filter = 90

    result_select_fasta = open("{}_selected_orfs_annotated_pts.fasta".format(output_path), "w")
    result_select_fasta_two =  open("{}_selected_orfs_annotated_pts_low_ident.fasta".format(output_path), "w")
    result_select_fasta_partial = open("{}_selected_orfs_annotated_pts_partial.fasta".format(output_path), "w")

    dict_blast_hit = {}

    for r in blastp_tab:
        report = r.strip()
        hit_fields = report.split('\t')
        transcript = hit_fields[0]
        toxin = hit_fields[3]
        ident= hit_fields[6]
        qstart = hit_fields[10]
        qend = hit_fields[11]
        sstart = hit_fields[12]
        send = hit_fields[13]
        coord_align = qstart + ":" + qend + "-" + hit_fields[1]
        coord_subject = sstart + ":" + send + "-" + hit_fields[4]
        cov = (int(hit_fields[7]) / int(hit_fields[4])) *100
        info = toxin + "|" + str(cov) + "|" + ident + "|" + coord_align + "|" + coord_subject

        if transcript not in dict_blast_hit.keys():
            dict_blast_hit[transcript] = info

    for k,v in dict_blast_hit.items():
        toxin, cov, ident, start_end_alig , start_end_subj = v.split("|")

        if (float(cov) >= float(cov_filter)) and (float(ident) >= 50.0):
            input_fasta = open(orfs_fasta, "r")
            fasta_sequences = SeqIO.parse(input_fasta,'fasta')

            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                len_seq = len(sequence)
                seq_id = name.split(" ")[0]

                if k == seq_id:
                    result_select_fasta.write(">{}|{}|{}|{}\n{}\n".format(k,toxin,start_end_alig,start_end_subj,sequence))

        elif (float(cov) >= float(cov_filter)) and (float(ident) <= 50.0):
            input_fasta = open(orfs_fasta, "r")
            fasta_sequences = SeqIO.parse(input_fasta,'fasta')

            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                len_seq = len(sequence)
                seq_id = name.split(" ")[0]

                if k == seq_id:
                    result_select_fasta_two.write(">{}|{}|{}|{}\n{}\n".format(k,toxin,start_end_alig,start_end_subj,sequence))

        elif (float(cov) >= 50.0) and (float(cov) <= 89.99) and (float(ident) >= 50.0):
            input_fasta = open(orfs_fasta, "r")
            fasta_sequences = SeqIO.parse(input_fasta,'fasta')

            for fasta in fasta_sequences:
                name, sequence = fasta.id, str(fasta.seq)
                len_seq = len(sequence)
                seq_id = name.split(" ")[0]

                if k == seq_id:
                    result_select_fasta_partial.write(">{}|{}|{}|{}\n{}\n".format(k,toxin,start_end_alig,start_end_subj,sequence))


    blastp_tab.close()
    input_fasta.close()
    result_select_fasta.close()
    result_select_fasta_two.close()
    result_select_fasta_partial.close()

  def get_ids_nts_orfs(self, out_dir_name, sample_name):
    os.system("grep '>' {}/{}_selected_orfs_annotated_pts.fasta | cut -f 1,2 -d '|' | sed 's/>//g' > {}/ids_get_nts_orfs"\
          .format(out_dir_name, sample_name, out_dir_name))

  def format_nts_orfs(self, out_dir_name, sample_name):
    os.system("sed -i 's/\s/|/g' {}/{}_selected_contigs_from_tblastn.fasta.transdecoder_dir/longest_orfs.cds"\
          .format(out_dir_name, sample_name))

  def get_nts_orfs(self, ids_nts_path, orfs_nts_path, out_dir_name, sample_name):
    list_to_get_nts = open(ids_nts_path, "r")
    orfs_nts_path = orfs_nts_path
    outfile_nts = open("{}/{}_selected_orfs_annotated_nts.fasta".format(out_dir_name, sample_name), "w")

    mydict = {}

    for i in list_to_get_nts:
        list_item = i.strip()
        orf_id , annotation = list_item.split("|")
        input_fasta = open(orfs_nts_path, "r")
        fasta_sequences = SeqIO.parse(input_fasta,'fasta')

        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            seq_id = name.split("|")[0]
            coord = name.split(":")[-1]
            header = orf_id + "|" + coord + "|" + annotation
            if orf_id == seq_id:
                mydict[header]= sequence
        input_fasta.close()

    for k,v in mydict.items():
        outfile_nts.write(">{}\n{}\n".format(k,v))

    list_to_get_nts.close()
    outfile_nts.close()

  def remove_redundancy(self, toxin_nts_path, toxin_pts_path):
    nts_file = toxin_nts_path
    pts_file = toxin_pts_path
    list_files = [nts_file,pts_file]

    for f in list_files:
        dedup_records = defaultdict(list)

        for record in SeqIO.parse(f, "fasta"):
            # Use the sequence as the key and then have a list of id's as the value
            dedup_records[str(record.seq)].append(record.id)
            out_name = f.replace(".fasta", "")
        with open(out_name + '_uniq.fasta', 'w') as output:
            for seq, ids in dedup_records.items():
                # Join the ids and write them out as the fasta
                output.write(">{}\n".format(ids[0]))
                output.write(seq + "\n")

  def get_ids_whole_contigs(self, out_dir_name, sample_name):
    os.system("grep '>' {}/{}_selected_orfs_annotated_nts_uniq.fasta | cut -f 1,2,3 -d '|' | sed 's/>//g' > {}/ids_get_whole_contigs"\
          .format(out_dir_name, sample_name, out_dir_name))

  def get_whole_contigs_from_selected_orfs(self, out_dir_name, sample_name):
    list_to_get_contigs = open("{}/ids_get_whole_contigs".format(out_dir_name), "r")
    outfile_whole_contigs = open("{}/{}_selected_whole_contigs.fasta".format(out_dir_name, sample_name), "w")
    whole_contigs = "{}/{}_selected_contigs_from_tblastn.fasta".format(out_dir_name, sample_name)
    mydict = {}

    for i in list_to_get_contigs:
        list_item = i.strip()
        contig_id , coord, annotation = list_item.split("|")
        contig_id_clean = re.sub('\.p.*','', contig_id)
        header = '{}|{}|{}'.format(contig_id_clean,coord,annotation)
        input_fasta = open(whole_contigs, "r")
        fasta_sequences = SeqIO.parse(input_fasta,'fasta')

        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            seq_id = name.split(" ")[0]

            if contig_id_clean == seq_id:
                mydict[header]= sequence
        input_fasta.close()

    for k,v in mydict.items():
        outfile_whole_contigs.write(">{}\n{}\n".format(k,v))

    list_to_get_contigs.close()
    outfile_whole_contigs.close()

  def count_toxin_families(self, fasta_file_path, out_dir_name, sample_name):
    input_file = open(fasta_file_path, "r")
    fasta_sequences = SeqIO.parse(input_file,'fasta')
    result = open("{}/{}_total_toxin_families.txt".format(out_dir_name, sample_name), "w")
    mydict = {}
    all_toxins = ["3FTx", "5NUCL", "ACES", "BDEF", "BPP", "CNP", "CRISP", "CTL", "CVF", "CYS", "DIESTER",
                "DIPEP", "DIS", "FA5V", "FAXV", "HYAL", "IPLA2", "KUNZ", "KUWAP", "LAO", "NGF", "OHAN",
                "PLA2", "PLB", "SBPM", "SRTX", "SVLIPA", "SVMI", "SVMMP", "SVMP", "SVSP", "VEGF", "WAP"]

    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        toxin = name.split("|")[-1]
        family = toxin.split("_")[-1]
        if family not in mydict.keys():
            mydict[family] = 1
        else:
            mydict[family] += 1
    for t in all_toxins:
        if t in mydict.keys():
            continue
        else:
            mydict[t] = 0
    for k,v in sorted(mydict.items()):
        result.write("{}\t{}\n".format(k,v))
    input_file.close()
    result.close()
