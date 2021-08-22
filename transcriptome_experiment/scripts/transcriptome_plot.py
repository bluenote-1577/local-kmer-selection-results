import numpy as np
import sys
import glob
import matplotlib.pyplot as plt
from natsort import natsorted
import re
from collections import defaultdict
import scipy.stats
import pickle
import matplotlib
plt.style.use('science')

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return h

iterate_dirs = sys.argv[3:]
num_iterates = len(iterate_dirs)
mapq_cutoff = 1
transcriptome_ref = sys.argv[1]
densmeth_to_timelist = pickle.load(open(sys.argv[2], 'rb'))
print(densmeth_to_timelist)

plt.rcParams.update({'font.size': 14})
#densities = ['d4', 'd6']
densities = ['d3', 'd5', 'd7']
xran = range(11,26)

gene_to_loc = dict()
for line in open(transcriptome_ref):
    if line[0] == '>':
        line_chop = line[1:]
        sp = line_chop.split()
        gene_name = sp[0].split('.')[0] 
        loc = sp[2].split(':')
        try:
            chrom_loc = int(loc[2])
        except ValueError:
            continue
        start_pos_on_chrom = int(loc[3])
        end_pos_on_chrom = int(loc[4])
        gene_to_loc[gene_name] = [chrom_loc,start_pos_on_chrom,end_pos_on_chrom]

sucs_sync_dict = dict()
sucs_mini_dict = dict()
fail_mini_dict = dict()
fail_sync_dict = dict()

for density in densities:
    sucs_order_sync_total = []
    sucs_order_mini_total = []
    fail_mini_total = []
    fail_sync_total = []

    samples_point_sucs_sync = []
    samples_point_sucs_mini = []
    samples_point_fail_sync = []
    samples_point_fail_mini = []

    for i in xran:
        samples_point_sucs_sync.append([])
        samples_point_sucs_mini.append([])
        samples_point_fail_sync.append([])
        samples_point_fail_mini.append([])


    for direct in iterate_dirs:
        sucs_fchrom_fpose = dict()
        sam = glob.glob(direct+'/*')
        for samfile in sam:
            if density not in samfile:
                continue
            if 'k10' in samfile:
                continue
            suc = 0
            fail_chrom = 0
            fail_pos = 0
            not_in = 0
            for line in open(samfile):
                if line[0] == '@':
                    continue
                splitted = line.split('\t')
                name = splitted[0]
                flag = int(splitted[1])
                ref = splitted[2]
                aln_pos = int(splitted[3])
                mapq = int(splitted[4])

                if flag & 256 or flag & 2048 or flag & 4:
                    continue

                if mapq < mapq_cutoff:
                    continue


                name_truth = name.split('_')[0]
                if name_truth not in gene_to_loc:
                    #print('not in', name_truth)
                    not_in += 1
                    continue
                start_aln_pos_truth = int(name.split('_')[1])
                map_length_truth = int(name.split('_')[-2]) + int(name.split('_')[-3]) + int(name.split('_')[-1])
                name_ref_map = int(ref.split('.')[0][3:])
                identi = ref.split('.')[0][0:2]
                if identi != "NC":
        #            print(name)
                    continue
                
                if gene_to_loc[name_truth][0] != name_ref_map:
                    fail_chrom += 1
            #        print(name_ref_map,aln_pos, gene_to_loc[name_truth])
            #        print(name,ref)
                elif aln_pos  > gene_to_loc[name_truth][2] or aln_pos + map_length_truth < gene_to_loc[name_truth][1]:
                    #print(aln_pos,gene_to_loc[name_truth][1] + start_aln_pos_truth,aln_pos - (gene_to_loc[name_truth][1] + start_aln_pos_truth), name)
                    fail_pos +=1
                else:
                    suc += 1

            print(samfile,suc,fail_chrom,fail_pos,  not_in, 'suc fail_chrom, fail_pos, not_in')
            sucs_fchrom_fpose[samfile] = (suc, fail_chrom, fail_pos)
                    #print(aln_pos, gene_to_loc[name_truth][1] + start_aln_pos_truth)
                #s1 = int(splitted[-4].split(':')[-1])
                #s2 = int(splitted[-3].split(':')[-1])


        sucs_order_sync = []
        sucs_order_mini = []
        fchrom_sync = []
        fchrom_mini= []
        fpos_sync = []
        fpos_mini= []
        fail_mini = []
        fail_sync = []

        for (i,key) in enumerate(natsorted(sucs_fchrom_fpose)):
            if 'mini' in key:
                sucs_order_mini.append(sucs_fchrom_fpose[key][0])
                fchrom_mini.append(sucs_fchrom_fpose[key][1])
                fpos_mini.append(sucs_fchrom_fpose[key][2])
            else:
                sucs_order_sync.append(sucs_fchrom_fpose[key][0])
                fchrom_sync.append(sucs_fchrom_fpose[key][1])
                fpos_sync.append(sucs_fchrom_fpose[key][2])

        fail_sync = np.array(fchrom_sync) + np.array(fpos_sync)
        fail_mini = np.array(fchrom_mini) + np.array(fpos_mini)

        for i in range(len(fail_sync)):
            samples_point_sucs_mini[i].append(sucs_order_mini[i])
            samples_point_fail_mini[i].append(fail_mini[i])
            samples_point_sucs_sync[i].append(sucs_order_sync[i])
            samples_point_fail_sync[i].append(fail_sync[i])



        if len(sucs_order_sync_total) == 0:
            sucs_order_sync_total = np.zeros(len(sucs_order_sync))
            sucs_order_mini_total = np.zeros(len(sucs_order_sync))
            fail_mini_total = np.zeros(len(sucs_order_sync))
            fail_sync_total = np.zeros(len(sucs_order_sync))
        print(sucs_order_sync_total,sucs_order_sync)
        sucs_order_sync_total += np.array(sucs_order_sync)/ num_iterates
        sucs_order_mini_total += np.array(sucs_order_mini)/ num_iterates
        fail_mini_total += np.array(fail_mini)/ num_iterates
        fail_sync_total += np.array(fail_sync)/ num_iterates


    fail_mini_dict[density] = fail_mini_total
    fail_sync_dict[density] = fail_sync_total
    sucs_sync_dict[density] = sucs_order_sync_total
    sucs_mini_dict[density] = sucs_order_mini_total

    conf_intervals_suc_sync = [mean_confidence_interval(x) for x in samples_point_sucs_sync]
    conf_intervals_suc_mini = [mean_confidence_interval(x) for x in samples_point_sucs_mini]
    conf_intervals_fail_sync = [mean_confidence_interval(x) for x in samples_point_fail_sync]
    conf_intervals_fail_mini = [mean_confidence_interval(x) for x in samples_point_fail_mini]

    print(conf_intervals_suc_sync)
    d5_marker = '--'
    d3_marker = '-'
    d7_marker = ':'
    mini_colour = 'C3'
    sync_colour = 'C0'

    if density == 'd5' or density == 'd6':
        marker = d5_marker
    elif density == 'd7':
        marker = d7_marker
    else:
        marker = d3_marker 

    plt.subplot(221)
    plt.errorbar(xran,sucs_order_sync_total, ls = marker,color = sync_colour, yerr = conf_intervals_suc_sync, label='Open syncmer, d = 1/' + density[1])
    plt.errorbar(xran, sucs_order_mini_total,ls = marker, color = mini_colour, yerr = conf_intervals_suc_mini, label='Minimizer, d = 1/' + density[1])
    plt.xlabel("k")
    plt.ylabel("number of reads successfully mapped")
    plt.legend()
    plt.subplot(223)
    plt.errorbar(xran,fail_sync_total, ls = marker, color = sync_colour, yerr = conf_intervals_fail_sync)
    plt.errorbar(xran,fail_mini_total, ls = marker, color = mini_colour, yerr = conf_intervals_fail_mini)
    plt.xlabel("k")
    plt.ylabel("number of errors")
    plt.subplot(222)
    plt.plot(sucs_order_sync_total, fail_sync_total, marker, color = sync_colour)
    plt.plot(sucs_order_mini_total, fail_mini_total, marker, color = mini_colour)
    plt.xlabel("number of reads successfully mapped")
    plt.ylabel("number of errors")
    plt.subplot(224)
    value_d = density[1]
    print(value_d+'sync', densmeth_to_timelist[value_d+'sync'])
    plt.plot(xran,densmeth_to_timelist[value_d+'sync'],marker, color = sync_colour)
    plt.plot(xran,densmeth_to_timelist[value_d+'mini'],marker, color = mini_colour)
    plt.xlabel("k")
    plt.ylabel("CPU time (s)")
    plt.yscale("log")

    
plt.show()


#plt.subplot(232)
for density in ['d3','d5','d7']:
    d5_marker = '--'
    d3_marker = '-'
    d7_marker = ':'
    mini_colour = 'C3'
    sync_colour = 'C0'

    if density == 'd5' or density == 'd6':
        marker = d5_marker
    elif density == 'd7':
        marker = d7_marker
    else:
        marker = d3_marker 


    value_d = density[1]
    plt.plot( fail_sync_dict[density],  densmeth_to_timelist[value_d+'sync'], marker+'o', markersize = 4, color = sync_colour, label='Open syncmer, d = 1/' + density[1])
    plt.plot(  fail_mini_dict[density], densmeth_to_timelist[value_d+'mini'], marker+'o', markersize = 4, color = mini_colour,  label='Minimizer, d = 1/' + density[1])   
    plt.xlabel("number of errors")
    plt.ylabel("CPU time (s)")
plt.legend()
plt.show()

#plt.subplot(235)
#plt.plot(  sucs_order_sync_total, densmeth_to_timelist[value_d+'sync'], marker+'o', markersize = 4, color = sync_colour)
#plt.plot(  sucs_order_mini_total, densmeth_to_timelist[value_d+'mini'], marker+'o', markersize = 4, color = mini_colour)
#plt.xlabel("number of reads successfully mapped")
#plt.ylabel("time")
#plt.yscale('log')



