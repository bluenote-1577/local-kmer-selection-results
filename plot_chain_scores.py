import re
import sys
import matplotlib
#matplotlib.use('') 
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
plt.style.use('science')
matplotlib.rc('text', usetex = False)

print(sys.argv)
sam = sys.argv[1:]
scores_sam = []
for (i,sam_file) in enumerate(sam):
    scores_sam.append([])
    scores = [re.findall(r's1:i:(\d+)',line) for line in open(sam_file)]
    for score in scores:
        if len(score) > 0:
            int_score = int(score[0])
            scores_sam[i].append(int_score)
#plt.plot([1,2,3])
colour = ['C0','C3']
for score in scores_sam:
    print(np.mean(score));
labels = ["Open syncmer, (k,s,t) = (15,11,3). Mean score= %d"%(np.mean(scores_sam[0])), "Random minimizer, (k,w) = (15,9). Mean score = %d"%(np.mean(scores_sam[1]))]
plt.figure(figsize=(17,12),dpi=80)
plt.rcParams.update({'font.size': 18})
for i in range(len(scores_sam)):
    plt.hist(scores_sam[i],color=colour[i],bins=300,density=True,alpha=0.35,label=labels[i])
    plt.xlim((0,15000))
    #plt.ylim((0,0.0005))
plt.legend()
plt.title("Chaining scores for long-reads from E. coli W (bc1087) mapping to assembly E. coli W (bc1087)")
plt.xlabel("Chaining score")
plt.ylabel("Density")
plt.show()


