import pysam
import re
import os
import collections
import sys
import getopt
import numpy as np
def get_best40(bam_file,f_dir_out):
    print(bam_file)
    best_40bam=bam_file.replace(".bam","best40.bam")
    r1="samtools view -h -F 0x100 '"+bam_file+"' | awk -F \'\\t\' \'OFS=\"\\t\" {gsub(/\/[12]$/,\"\",$1);{print}}\' - |samtools view -h -Sb -| samtools rmdup - - | samtools view -h -b -F 4 -q 40 -  >'"+best_40bam+"' \n"
    r2="fastqc -o '"+f_dir_out+"' --extract -t 10 -q -d '"+f_dir_out+"' '"+best_40bam+"'\n"
    if(os.path.exists(best_40bam)):
        if(os.path.getsize(best_40bam)<10000000000):
            print(best_40bam)
            os.system(r1)
    else:
        os.system(r1+r2)

def fetch_jun(st, cigar):
    """ fetch exon regions defined by cigar. st must be zero based
	return list of tuple of (chrom,st, end)
	"""
    # match = re.compile(r'(\d+)(\D)')
    chrome_st = st
    jun_bound = []
    jun_read_pos = []
    read_st = 0
    chrome_e = chrome_st
    read_e = read_st
    for c, s in cigar:  # code and size
        if c == 0:  # match
            chrome_e += s
            read_e += s
        elif c == 1:  # insertion to ref
            read_e += s
        elif c == 2:  # deletion to ref
            chrome_e += s
        elif c == 3:  # gap or intron
            jun_bound.append([chrome_st, chrome_e])
            jun_read_pos.append([read_st, read_e])
            chrome_e += s
            chrome_st = chrome_e
            read_st = read_e
        else:
            continue
    jun_bound.append([chrome_st, chrome_e])
    jun_read_pos.append([read_st, read_e])
    return [jun_bound, jun_read_pos]


def filter_current_junction(bam_file, jun_file, f_result):  # pysam is 0 base coordination
    # ############  store the junction information in a dist
    f_jun = open(jun_file, "r")
    lines = f_jun.readlines()
    f_jun.close()
    all_jun_dist = collections.defaultdict(dict)
    for line in lines:
        temp = line.split("\t")
        all_jun_dist[temp[0]][",".join(temp[1:3])] = [temp[3], 0]
    # ################## read bamfile
    print("read bam")
    bam = os.path.split(str(bam_file))
    bam_name = bam[-1].split(".")
    print(bam_name[0])
    f_w_n = open(f_result + bam_name[0] + "_all_new_junction_read.bed", "w")
    f_w_c = open(f_result + bam_name[0] + "_all_current_read_count.bed", "w")
    f_w_s = open(f_result + bam_name[0] + "_statistical_num.txt", "w")
    f_w_in = open(f_result + bam_name[0] + "_good_reads.txt", "w")
    f_w_nc = open(f_result + bam_name[0] + "_all_new_junction_count.txt", "w")
    bam_reader = pysam.Samfile(bam_file, "rb")
    bam_write = pysam.Samfile(f_result + bam_name[0] + "_all_abnormal_junction_reads.bam", "wb", template=bam_reader)
    good_reads_file=f_result + bam_name[0] + "_good_reads.txt"
    chr_array = bam_reader.references
    reads_all_count = 0
    jun_read_all_count = 0
    new_jun_read_count = 0  # number of new junction reads count
    new_junction_count = 0  # number of new junction that has overlapped reads
    all_junction = 0
    if chr_array[0][0:3] == "chr":
        flag_chr = 3
    else:
        flag_chr = 0
    for chr_i in chr_array:
        chr_j = chr_i[flag_chr:]
        if not re.match("^[0-9xyXY]+$", chr_j):
            continue
        print(chr_j)
        chr_jun = all_jun_dist[chr_j]
        chr_reads = bam_reader.fetch(chr_i)
        all_read_junction = {}
        for read_seq in chr_reads:
            reads_all_count += 1
            cigar = read_seq.cigarstring
            if re.findall("N", cigar):  # d="4M11250N5M1D39M44N99M"
                jun_read_all_count += 1
                jun_all = fetch_jun(read_seq.pos, read_seq.cigar)
                i = 0
                chr_pos = jun_all[0]
                read_pos = jun_all[1]
                flag = "not new"
                flag_good="F"
                while i < len(chr_pos) - 1:
                    new_junction = str(int(chr_pos[i][1]) - 1) + "," + str(chr_pos[i + 1][0] + 1)
                    if new_junction in chr_jun:
                        chr_jun[new_junction][1] += 1
                    else:
                        junction_quality = ord(read_seq.qual[read_pos[i][1] - 2]) + ord(
                            read_seq.qual[read_pos[i][1] - 1]) + ord(read_seq.qual[read_pos[i + 1][0] - 1]) + ord(
                            read_seq.qual[read_pos[i + 1][0]])
                        average_junction_quality = (junction_quality - 33 * 4 - 30 * 4) / 4.0
                        if chr_j + "," + new_junction in all_read_junction:
                            all_read_junction[chr_j + "," + new_junction][0] += 1
                            all_read_junction[chr_j + "," + new_junction][1].append(chr_pos[i][0] + 1)
                            all_read_junction[chr_j + "," + new_junction][2].append(chr_pos[i + 1][1])
                            all_read_junction[chr_j + "," + new_junction][3].append(average_junction_quality)
                            all_read_junction[chr_j + "," + new_junction][4].append(read_seq.isize)
                            all_read_junction[chr_j + "," + new_junction][5].append(read_seq.tags[-1][1])
                            all_read_junction[chr_j + "," + new_junction][6].append(read_seq.flag)
                            all_read_junction[chr_j + "," + new_junction][7].append(read_seq.qname)
                        else:
                            all_read_junction[chr_j + "," + new_junction] = [1, [chr_pos[i][0] + 1],
                                                                             [chr_pos[i + 1][1]],
                                                                             [average_junction_quality],
                                                                             [read_seq.isize],
                                                                             [read_seq.tags[-1][1]], [read_seq.flag],[read_seq.qname]]
                        print >> f_w_n, "%s\t%d\t%d\t%d\t%d\t%s,%s,%s,%s\t%s\t%s\t%s\t%s" % (
                            chr_j, chr_pos[i][0] + 1, chr_pos[i][1], chr_pos[i + 1][0] + 1, chr_pos[i + 1][1],
                            read_seq.qname,
                            read_seq.cigarstring, read_seq.mapq, read_seq.flag,
                            read_seq.qual[read_pos[i][0]:read_pos[i][1]],
                            read_seq.qual[read_pos[i + 1][0]:read_pos[i + 1][1]], read_seq.isize, read_seq.tags[-1][1]
                        )
                        flag = "new"
                        if average_junction_quality >0 and chr_pos[i][1] - chr_pos[i][0]>15 and chr_pos[i+1][1] - chr_pos[i+1][0]>15:
                            flag_good="T"
                    if flag == "new" and flag_good=="T":
                        f_w_in.write(chr_j+"\t"+read_seq.qname+"\t"+str(chr_pos[i][1])+";"+str(chr_pos[i + 1][0] + 1)+"\n")
                        new_jun_read_count += 1
                        bam_write.write(read_seq)
                    i += 1
        for a_chr_jun in chr_jun:
            if chr_jun[a_chr_jun][1] > 0:
                all_junction += 1
            print >> f_w_c, "%s\t%s\t%s\t%d" % (
                chr_j, a_chr_jun.replace(",", "\t"), chr_jun[a_chr_jun][0], chr_jun[a_chr_jun][1])
        for a_r_jun in all_read_junction:
            temp_a = a_r_jun.split(",")
            t_array = all_read_junction[a_r_jun]
            temp_left = ";".join(map(str, t_array[1]))
            temp_right = ";".join(map(str, t_array[2]))
            temp_score = ";".join(map(str, t_array[3]))
            print >> f_w_nc, "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s" % (
                temp_a[0], min(t_array[1]), int(temp_a[1]) + 1, temp_a[2], max(t_array[2]),
                t_array[0], str(max(t_array[3])), temp_left, temp_right, temp_score)
            new_junction_count += 1
    print >> f_w_s, "bam:%s\tall_read:%d\tall_junction_read:%d\tnew_junction_read:%d\tall_current_junction_with_read_count:%d\tnew_junction_with_read_count:%d" % (
        bam_name[0], reads_all_count, jun_read_all_count, new_jun_read_count, all_junction, new_junction_count)
    f_w_n.close()
    f_w_c.close()
    f_w_s.close()
    f_w_nc.close()
    bam_reader.close()
    bam_write.close()
    f_w_in.close()
    os.system("samtools index " + f_result + bam_name[0] + "_all_abnormal_junction_reads.bam")
    return good_reads_file

def get_mate_junction(bam_file,  f_result):  # pysam is 0 base coordination
    # ############  store the junction information in a dist

    # ################## read bamfile
    bam = os.path.split(str(bam_file))
    bam_name = bam[-1].split(".")
    jun_file=f_result + bam_name[0] + "_good_reads.txt"
    f_jun = open(jun_file, "r")
    lines = f_jun.readlines()
    f_jun.close()
    print(len(lines))
    all_jun_dist = collections.defaultdict(dict)
    for line in lines:
        temp = line.strip().split("\t")
        all_jun_dist[temp[0]][temp[1]] = [temp[2]]
    print(bam_name[0])
    f_w_s = open(f_result + bam_name[0] + "_non_jun_reads_isize_num.txt", "w")
    f_w_g = open(f_result + bam_name[0] + "_good_reads_isize_num.txt", "w")
    f_w_c = open(f_result + bam_name[0] + "_good_reads_count.txt", "w")
    f_w_in = open(f_result + bam_name[0] + "_good_reads_quantile.txt", "w")
    bam_reader = pysam.Samfile(bam_file, "rb")
    bam_write = pysam.Samfile(f_result + bam_name[0] + "_all_new_junction_with_mate.bam", "wb", template=bam_reader)
    chr_array = bam_reader.references
    flag = [83, 99, 147, 163]
    insert_size_array=[]
    if chr_array[0][0:3] == "chr":
        flag_chr = 3
    else:
        flag_chr = 0
    for chr_i in chr_array:
        chr_j = chr_i[flag_chr:]
        if not re.match("^[0-9xyXY]+$", chr_j):
            continue
        print(chr_j)
        chr_jun = all_jun_dist[chr_j]
        print(len(chr_jun))
        chr_reads = bam_reader.fetch(chr_i)
        print(chr_reads.next().qname)
        jun_count={}
        i = 0
        for read_seq in chr_reads:
            if read_seq.flag not in flag:
                continue
            r_name = read_seq.qname.split("/")[0]
            if r_name in chr_jun:
                chr_jun[r_name].append(read_seq)
            if read_seq.isize > 0 and not re.findall("N", read_seq.cigarstring) and i < 500:
                f_w_s.write(r_name + "\t" + str(read_seq.mpos - (read_seq.pos + read_seq.cigar[0][1])) + "\n")
                i = i + 1
        for r in chr_jun:
            if len(chr_jun[r]) > 2:
                r1 = chr_jun[r][1]
                r2 = chr_jun[r][2]
                if r1.pos < r2.pos:
                    first_r = r1
                    second_r = r2
                else:
                    first_r = r2
                    second_r = r1
                if re.findall("N", first_r.cigarstring):
                    jun_all = fetch_jun(first_r.pos, first_r.cigar)
                    chr_pos = jun_all[0]
                    insert_size =  second_r.pos - chr_pos[-1][1]
                else:
                    insert_size = -(first_r.pos + first_r.cigar[0][1]) + second_r.pos
                if insert_size>35:
                    bam_write.write(first_r)
                    bam_write.write(second_r)
                    if chr_jun[r][0] in jun_count:
                        jun_count[chr_jun[r][0]].append(r)
                    else:
                        jun_count[chr_jun[r][0]]=[r]
                f_w_g.write(r + "\t" + str(insert_size) + "\n")
                insert_size_array.append(insert_size)
            if len(chr_jun[r]) > 3:
                print(chr_jun[r])
        for jun in jun_count:
            f_w_c.write(chr_j+";"+jun+"\t"+str(len(jun_count[jun]))+"\t"+";".join(jun_count[jun])+"\n")
    insert_np=np.array(insert_size_array)
    p_t_10 = np.percentile(insert_np, 10)
    p_t_50 = np.percentile(insert_np, 50)
    p_t_90 = np.percentile(insert_np, 90)
    f_w_in.write(bam_name[0]+"\n"+str(p_t_10)+"\n"+str(p_t_50)+"\n"+str(p_t_90)+"\n")
    f_w_s.close()
    f_w_c.close()
    bam_reader.close()
    bam_write.close()
    f_w_g.close()
    f_w_in.close()
    os.system("samtools sort " + f_result + bam_name[0] + "_all_new_junction_with_mate.bam   "+ f_result + bam_name[0] + "_all_new_junction_with_mate_sort")
    os.system("samtools index " + f_result + bam_name[0] + "_all_new_junction_with_mate_sort.bam")


def run_all():
    opts, args = getopt.getopt(sys.argv[1:], "hb:o:r:")
    f_out = ""
    current_jun = ""
    filename = ""
    print(opts)
    for op, value in opts:
        if op == "-b":
            filename = value
        elif op == "-o":
            f_out = value
        elif op == "-r":
            current_jun = value
        elif op == "-h":
            print("please input the parameters\n")
            sys.exit()
    f_result = f_out + '/step1_filter_current_junction_read/'
    if not os.path.exists(f_result):
        os.mkdir(f_result)
    best_40=filename.replace(".bam","best40.bam")
    if not os.path.isfile(best_40):
        get_best40(filename.strip(),f_result)
        os.system("samtools index " + best_40)
    jun_file=filter_current_junction(best_40, current_jun, f_result)
    get_mate_junction(best_40,  f_result)


if __name__ == '__main__':
    run_all()
