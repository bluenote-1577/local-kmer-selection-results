#Look at reads to see if they align using one method, but not the other. 
import sys
import re
from collections import defaultdict
sam = sys.argv[1:]
reads_map_methods = []
for i in range(len(sam)):
    reads_map_methods.append(defaultdict(list))
for (i,sam_file) in enumerate(sam):
    reads_map = reads_map_methods[i]
    for line in open(sam_file):
        if line[0] != 'm':
            continue
        splitted = line.split('\t')
        name = splitted[0]
        flag = int(splitted[1])
        mapq = int(splitted[4])
        reads_map[splitted[0]].append([flag,mapq])

unmapped_reads = set()
for read in reads_map_methods[0].keys():
    for i in range(len(sam)):
        if reads_map_methods[i][read][0][0] == 4:
            unmapped_reads.add(read)

#print(unmapped_reads)

sync_map = 0
min_map = 0
avg_mapq_sync = 0
avg_mapq_min = 0
for read in unmapped_reads:
    if reads_map_methods[0][read][0][0] != 4:
 #       print(reads_map_methods[0][read],reads_map_methods[1][read],"SYNC MAP")
        sync_map +=1
        max_mapq = max([x[1] for x in reads_map_methods[0][read]])
        avg_mapq_sync += max_mapq
    elif reads_map_methods[1][read][0][0] != 4:
 #       print(reads_map_methods[0][read],reads_map_methods[1][read],"MIN MAP")
        max_mapq = max([x[1] for x in reads_map_methods[1][read]])
        min_map +=1
        avg_mapq_min += max_mapq

print(min_map,sync_map,len(reads_map_methods[0].keys()), len(unmapped_reads),"size (M and not OS), size (OS and not M), num reads, num unmapped reads (R)")
print(avg_mapq_min/min_map, avg_mapq_sync/sync_map,"avg mapq syncmer, avg mapq minimizer")
