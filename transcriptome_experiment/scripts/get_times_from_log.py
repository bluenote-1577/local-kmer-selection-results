import re
import numpy as np
import pickle
from natsort import natsorted
from collections import defaultdict
big_log = '/home/jshaw/scratch/2021_theory_kmers_human_experiment/all_runs.log'

mean_densmeth_to_timelist_sorted = defaultdict(list)
iters = ['i1','i2','i3','i4','i5','i6','i7','i8','i9']
#iters = ['i2']
for it in iters:
    densmeth_to_timelist = defaultdict(dict)
    densmeth_to_timelist_sorted = defaultdict(list)

    f = open(big_log)
    for line in f:
        if "CMD" in line:
            if it in line:
                match = re.findall('k ([0-9]+)',line)
                value_k = int(match[0])
                if value_k < 11:
                    continue
                match = re.findall('d([0-9])',line)
                value_d = match[0]
                if 'syncs' in line:
                    densmeth_to_timelist[value_d + 'sync'][value_k] = int(re.findall('CPU: ([0-9]+)',next(f))[0])
                else:
                    densmeth_to_timelist[value_d + 'mini'][value_k] = int(re.findall('CPU: ([0-9]+)',next(f))[0])

    for densmeth in densmeth_to_timelist:
        for k in natsorted(densmeth_to_timelist[densmeth]):
            densmeth_to_timelist_sorted[densmeth].append(densmeth_to_timelist[densmeth][k])


    #print(densmeth_to_timelist_sorted)

    for key in densmeth_to_timelist_sorted:
        if key in mean_densmeth_to_timelist_sorted:
            print(mean_densmeth_to_timelist_sorted[key],densmeth_to_timelist_sorted[key])
            mean_densmeth_to_timelist_sorted[key] += np.array(densmeth_to_timelist_sorted[key])/len(iters)
        else:
            mean_densmeth_to_timelist_sorted[key] = np.array(densmeth_to_timelist_sorted[key])/len(iters)
    f.close()

print(mean_densmeth_to_timelist_sorted)
pickle.dump(mean_densmeth_to_timelist_sorted,open('dump.pkl','wb'))
