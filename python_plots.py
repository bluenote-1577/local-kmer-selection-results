import numpy as np
import scipy.special
import matplotlib.pyplot as plt
import matplotlib
plt.style.use('science')

matplotlib.rc('text', usetex = True)
run_prob_k_theta = False
syncmer_prob_vect_plots = False
prob_vect_compare_plots= False
cons_compare_plots = True

##Empirical results. Stored in emp_** text file, copy and pasted here. 
os_cons_over_theta = [
    0.8971449400000001,
    0.7854196600000003,
    0.6730598399999999,
    0.5673365799999999,
    0.47264378000000007,
    0.3891510400000001,
    0.31717677999999994,
    0.25556848000000004,
    0.20518112,
    0.16312411999999996,
    0.12895374000000004,
    0.10144844000000001,
    0.07910765999999998,
    0.06188425999999997,
    0.04722440000000001,
]
cept_dens_mean = 0.12158664365227469
min_cons_over_theta = [
    0.8745431800000002,
    0.7478939399999998,
    0.6276925999999999,
    0.5190208199999999,
    0.42436008000000003,
    0.3438545400000001,
    0.2754275400000001,
    0.2188662799999999,
    0.17304701999999994,
    0.13553726,
    0.10580112,
    0.08198255999999997,
    0.06346327999999998,
    0.04886084,
    0.036885259999999996,
]
min_cons_std_over_theta = [
    0.0016590120729463663,
    0.002782092773144671,
    0.0031788741352366846,
    0.0034787952877600365,
    0.0034442501219606935,
    0.0029792927078063185,
    0.003524908146722162,
    0.0028333787521493334,
    0.0030343628555250103,
    0.0025997370442785066,
    0.002241825200116813,
    0.0019389040970026702,
    0.0018450604221033312,
    0.0015560501713497467,
    0.0015096093723111607,
]
cept_cons_over_theta = [
    0.8830978999999999,
    0.7625268800000002,
    0.6456358799999998,
    0.53839866,
    0.4440552000000001,
    0.36223104,
    0.2925774,
    0.23418798000000007,
    0.18646187999999994,
    0.14716978000000003,
    0.11541270000000005,
    0.09030569999999997,
    0.07006085999999999,
    0.05441830000000001,
    0.041392620000000005,
]
cept_cons_std_over_theta = [
    0.0015771310963104794,
    0.002597297170428354,
    0.0029795925939394134,
    0.0036435304008917333,
    0.003349775665638076,
    0.0032077670698075047,
    0.0034527425901650106,
    0.0031969674582155026,
    0.003056551166913389,
    0.0026960895122717676,
    0.002414328614038622,
    0.001992019682167175,
    0.0017796324953725305,
    0.0015554721748937484,
    0.0015225512267177744,
]
cust_cons_over_theta = [
    0.8986747000000003,
    0.78735836,
    0.6754631600000002,
    0.56987722,
    0.47472906000000015,
    0.3910401199999999,
    0.3186199800000001,
    0.25715703999999995,
    0.20618859999999994,
    0.16427853999999997,
    0.12966049999999998,
    0.10223598000000003,
    0.07950638000000002,
    0.06212232000000004,
    0.04761830000000001,
]
cust_cons_std_over_theta = [
    0.0014183829283071325,
    0.0024669596448448246,
    0.002936060566628194,
    0.0035990048219871696,
    0.0033228796013169805,
    0.0033283469227200152,
    0.0035784671453999786,
    0.0030683130426468488,
    0.0032063099089643904,
    0.0027543526597263107,
    0.0024921647622138062,
    0.0021110612295117057,
    0.0018764976337029922,
    0.0018211582000185898,
    0.0016942755521061638,
]

#The scheme prob. vector for the "abbb..." words scheme from Frith(2020).
#n : number of "b"s
#k : k-mer k
def words_prob(n,k):
    vec_prob = []
    for a in range(k):
        prob_sum = 0
        for i in range(1, a+2):
            prob_sum += (-1)**(i+1) * 3**(n*i)/(4 * 4**n)**i * binom((a+1) - n*(i-1),i)
        vec_prob.append(prob_sum)

    return(vec_prob)

#closed ysncmer scheme
def closed_sync(k,s):
    prob_vec = []
    for i in range(1,k+1):
        alpha = k - s + i
        prob = 2 * i * np.math.factorial(alpha - 1)/np.math.factorial(alpha)
        if prob > 1:
            prob_vec.append(1)
        else:
            prob_vec.append(prob)
    return prob_vec

#The scheme prob. vector for the open-syncer with offset scheme from Edgar(2021).
#k : k-mer k
#s : s-mer s
#t : offset t
def union_os_prob_exact_middle(k,s,t):
    #L is the number of permutations which satisfy
    L_vec = [0,np.math.factorial(k-s)]
    prob_vec = [1/(k-s+1)]
    for i in range(2,k+1):
        alpha = k-s+i
        recur_term = 0
        for l in range(1,t):
            if i-l < 0:
                continue
            recur_term += L_vec[i-l] * np.math.factorial(alpha-1)/np.math.factorial(alpha-l)
        for l in range(1,k-s-t+2):
            if i-l < 0:
                continue
            recur_term += L_vec[i-l] * np.math.factorial(alpha-1)/np.math.factorial(alpha-l)

        L_vec.append(i*np.math.factorial(alpha-1) + recur_term)
        prob_vec.append(L_vec[-1]/np.math.factorial(alpha))
    
    #print(prob_vec)
    #print(L_vec)
    return prob_vec

def binom(n,k):
    if n < 0:
        return 0
    if n >= k:
        return int(scipy.special.binom(n,k))
    else:
        return 0

def profile(a_p1,theta,k,count = False):
    if a_p1 == k:
        return (1 - theta)**(2*k-1)

    a = a_p1 - 1
    profile_prob = 0

    for b in range(0,k-a-1):
        Tk = 2*binom(k-2-a,b) + (k-a-2)*binom(k-3-a,b)
        profile_prob += Tk * (1-theta)**(k+a+b)*(theta)**(k-a-b-1)
    if count:
        #if theta = 1/2, this returns a count instead of a probability
        return profile_prob*(2**(2*k-1))
    else:
        return profile_prob
    
def profile_vec(theta,k):
    profile_prob_vec = []
    for a_p1 in range(1,k+1):
        profile_prob_vec.append(profile(a_p1,theta,k))

    return(profile_prob_vec)

def falling_fac(a,b):
    term = 1
    for i in range(b,a):
        term *= (i+1)
    return term

#ideal minimizer scheme
def minimizer_prob_exact(w,k):
    return_vec = []
    for a in range(0,w):
        ind = w
        M = np.zeros((2*w-1+a,2*w-1))
        #M is the number of perms where M(A,B,C) is the defn given in the paper.

        #for p in range(1,(w+1)//2+1):
        for p in range(1,2*w):
            M[w-1,p-1] = (1+a)*np.math.factorial(w-1)
        for n in range(w+1,2*w+a):
            for p in range(1,2*w):
                M_term = (a+1)*np.math.factorial(n-1)
                l1 = p - 1
                r1 = n - (p + a)
                recur_term_2 = 0
                recur_term_1 = 0

                for i in range(1,l1+1):
                    recur_term_1 += M[n-i-1,l1 - i] * falling_fac(n-1,n-i)
                for i in range(1,r1+1):
                    #print(r1,l1,n,i,'r1 l1 n i',w,a)
                    recur_term_2 += M[n-i-1,r1 - i] * falling_fac(n-1,n-i)

                M[n-1,p-1] = M_term + recur_term_1 + recur_term_2

        np.set_printoptions(precision=2)
        return_vec.append(M[-1,w-1]/np.math.factorial(2*w-1+a))
    for i in range(k-w):
        return_vec.append(1.0)
    return return_vec


def union_bound(k,p):
    ret_vec = []
    for i in range(k):
        if 1/p*(i+1) > 1:
            ret_vec.append(1)
        else:
            ret_vec.append((i+1)/p)

    return ret_vec



if run_prob_k_theta:
    fig, axs = plt.subplots(3, 3)
    thetas = [0.01,0.05,0.10]
    ks = [10,15,20]

    for (i,theta) in enumerate(thetas):
        for (j,k) in enumerate(ks):
            profile_vector = profile_vec(theta,k)
            fig.suptitle(r'$\Pr(\alpha)$ for various values of $(\theta,k)$',fontsize=16)
            axs[i,j].plot(range(1,k+1),profile_vector,'o-')
            title = "(%s,%d)" % (theta,k)
            title= r'$(\theta,k)$ = ' + title
            axs[i,j].set_title(title)

    for ax in axs.flat:
        ax.set(xlabel=r'$\alpha$', ylabel='Probability')
    for ax in axs.flat:
        ax.label_outer()
    plt.show()

if syncmer_prob_vect_plots:
    k = 17
    s = 10
    cs_s = 2*s - k - 1
    ts = range(1,(k-s+2)//2+1)
    plt.figure(figsize=(13,11),dpi=80)
    plt.rcParams.update({'font.size': 18})

    for t in ts:
        os_vect = union_os_prob_exact_middle(k,s,t)
        label = r'Open syncmer with t = %s' %(t)
        plt.plot(range(1,k+1),os_vect,'-o',label = label)
    cs_vect = closed_sync(k,cs_s)
    plt.plot(range(1,k+1),cs_vect,'-o',label = 'Closed syncmer')
    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$\Pr(f,\alpha)$')
    plt.legend()
    plt.title(r'$\Pr(f)$ for syncmer methods with $(k,s) = (%s,%s)$'%(k,s))
    plt.show()

#Parameters for figure with k = 17 and comparing Pr(f) in the paper
#k = 17
#w = 13
#p = (w+1)//2
#s = k-p+1
#cs_s = 2*s - k - 1
#t = (k-s)//2 +1
#n = 2

#Parameters for k = 24 and comparing conservation in the paper
k = 24
w = 15
p = (w+1)//2
s = k-p + 1
cs_s = 2*s - k - 1
t = (k-s)//2 +1
n = 2
samp_size =100 


if prob_vect_compare_plots:
    
    os_exact = union_os_prob_exact_middle(k,s,t)
    words_exact = words_prob(n,k)
    mini_exact = minimizer_prob_exact(w,k)
    closed_sync = closed_sync(k,cs_s)
    ub = union_bound(k,p)

    plt.figure(figsize=(13,11),dpi=80)
    plt.rcParams.update({'font.size': 18})

    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$\Pr(f,\alpha)$')

    #colours for:
    #os:r
    #words:g
    #mini:b
    #cs:m
    #ub:k
    plt.title(r'Probability vectors for $k = %s$'%(k))
    plt.plot(range(1,k+1),os_exact,'-o',label = r'Open syncmer with optimal $t, d = 1/7$')
    plt.plot(range(1,k+1),words_exact,'-o',label= r'$(a,b,n)$-words method, $d = 9/64 = 1/7.11...$')
    plt.plot(range(1,k+1),mini_exact,'-o',label='Random minimizer, $d = 1/7$')
    plt.plot(range(1,k+1),closed_sync,'-o',label='Closed syncmer, $d = 1/7$')
    plt.plot(range(1,k+1),ub,'-o',label='Union bound, $d = 1/7$')
    plt.legend()
    plt.show()

if cons_compare_plots:

    thetas = np.linspace(0.01,0.15,15);
    mini_cons = []
    os_cons = []
    os_cons_emp = []
    words_cons = []
    words3_cons = []
    ub_cons = []
    cs_cons = []
    mini_cons_emp = []
    mini_cons_emp_er = []
    cept_cons_emp = []
    cept_cons_emp_er = []
    cust_cons_emp = []
    cust_cons_emp_er = []

    os_exact = union_os_prob_exact_middle(k,s,t)
    words_exact = words_prob(n,k)
    words3_exact = words_prob(n+1,k)
    mini_exact = minimizer_prob_exact(w,k)
    cs_exact = closed_sync(k,cs_s)
    ub = union_bound(k,p)

    plt.figure(figsize=(17,12),dpi=80)
    plt.rcParams.update({'font.size': 18})

    for (i,theta) in enumerate(thetas):
        profile_vector = profile_vec(theta,k)
        ub_cons.append(np.dot(ub,profile_vector))
        mini_cons.append(np.dot(mini_exact,profile_vector)/ub_cons[-1])
        os_cons.append(np.dot(os_exact,profile_vector)/ub_cons[-1])
        cs_cons.append(np.dot(cs_exact,profile_vector)/ub_cons[-1])
        words_cons.append(np.dot(words_exact,profile_vector)/ub_cons[-1])
        words3_cons.append(np.dot(words3_exact,profile_vector)/ub_cons[-1])
        mini_cons_emp.append(min_cons_over_theta[i]/ub_cons[-1])
        mini_cons_emp_er.append(min_cons_std_over_theta[i]/ub_cons[-1])
        cept_cons_emp.append(cept_cons_over_theta[i]/ub_cons[-1])
        cept_cons_emp_er.append(cept_cons_std_over_theta[i]/ub_cons[-1])
        cust_cons_emp.append(cust_cons_over_theta[i]/ub_cons[-1])
        cust_cons_emp_er.append(cust_cons_std_over_theta[i]/ub_cons[-1])
        os_cons_emp.append(os_cons_over_theta[i]/ub_cons[-1])

    mini_cons_emp_er = np.array(mini_cons_emp_er) * 1.92 / np.sqrt(samp_size);
    cept_cons_emp_er = np.array(cept_cons_emp_er) * 1.92 / np.sqrt(samp_size);
    cust_cons_emp_er = np.array(cust_cons_emp_er) * 1.92 / np.sqrt(samp_size);

    plt.plot(thetas,os_cons,'-o',label='Open syncmer with optimal $t$')
    plt.plot(thetas,words_cons,'-^',label='$(a,b,n)$-words, $d = 1/7.11$')
    plt.plot(thetas,words3_cons,'-s',c='C1',label='$(a,b,n)$-words, $d = 1/9.48$')
    plt.plot(thetas,mini_cons,'--o',c='C2',label=r'Random minimizer upper bound using $\Pr(f) \cdot \Pr(\alpha(\theta,k))$')
    plt.plot(thetas,cs_cons,'-o',c='C3',label='Closed syncmer')
    plt.errorbar(thetas,mini_cons_emp,fmt='-x',c='C2',yerr=mini_cons_emp_er,label=r'Empirical random minimizer')
    plt.errorbar(thetas,cept_cons_emp,fmt='-x',c='C4',yerr=cept_cons_emp_er,label=r'Empirical miniception $d = 0.1216, w = 13, k_0 = 11$')
    plt.errorbar(thetas,cust_cons_emp,fmt='-x',c='C6',yerr=cust_cons_emp_er,label=r'Empirical words on $W_8$')
    plt.plot(thetas,os_cons_emp,'-o',c='b',label='EMP OS')
    plt.xticks(thetas)
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\frac{\text{Cons}(f,t,\theta)}{UB(d) \cdot \Pr(\alpha(\theta,k))}$')
    plt.title(r'Computing $\frac{\text{Cons}(f,t,\theta)}{UB(d) \cdot \Pr(\alpha(\theta,k))}$ for $k = 24$ and $d \sim 1/8$')
    #plt.plot(ub_cons,label='ub cons')
    plt.legend()
    plt.show()

