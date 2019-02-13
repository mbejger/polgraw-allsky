#from itertools import combinations
from pygsl.combination import *
from scipy.special import binom 
import numpy as np 
import sys 


def FastAlarm(Cmax, noc, L, Nc, Nk): 

    ee = [x / Nc for x in Nk] 

    Nmax = min(noc, L) 

    c = [0 for _ in range(Nmax-Cmax+1)]
    pf = 0

    for i in range(Cmax, Nmax + 1): 
     
        cp = Combination(L, i)
        cq = Combination(L, L-i)
        cq.init_last()

        while True:
            p = 1
            for ind in cp.tolist():
                p *= ee[ind]
            #print(cp.tolist())

            q = 1
            for ind in cq.tolist():
                q *= 1 - ee[ind]
            #print(cq.tolist())

            c[i-Cmax] += p*q

            if cp.next() != 0:
                break
            if cq.prev() != 0:
                break


        pf += c[i-Cmax]  

    #r0 = 1 - np.power(1. - pf, Nc)
    #r1 = Nc*pf 
    #r2 = pf 

    #print pf, c 

    print(c)
    c.reverse()
    return c 


def main(): 

    Cmax = 2 
    noc = int(sys.argv[1])  

#    Nc = 956631135.0
#    Nk = [14707714, 28142054, 25801055, 24110051, 52838959, 28487013, 24250093, 19706750] 

    Nc = 305691837.0
    Nk = [5589016, 3514365, 4088149, 3859672, 4276426, 3399127, 3613982, 4091684, 4472328, 4057213, 5894944, 4420419, 4711823, 4224792, 3813295, 3672673, 3174463, 4399280, 2943254, 3888526, 3590689, 4212322, 3795185, 3645355, 3494263, 4597885]

    L = len(Nk) 
      
    C0 = FastAlarm(Cmax, noc, L, Nc, Nk) 
    C1 = FastAlarm(Cmax, noc, L, 2*Nc, Nk) 
    C2 = FastAlarm(Cmax, noc, L, 4*Nc, Nk) 
    C3 = FastAlarm(Cmax, noc, L, 8*Nc, Nk) 
    C4 = FastAlarm(Cmax, noc, L, 16*Nc, Nk) 

    print(C0)
    print(C1)
    print(C2)
    print(C3)
    print(C4)
    #print pf0, pf1, pf2, pf3, pf4
 
    pfe0 = np.cumsum(C0) 
    pfe1 = np.cumsum(C1) 
    pfe2 = np.cumsum(C2) 
    pfe3 = np.cumsum(C3) 
    pfe4 = np.cumsum(C4) 

    #print pfe0, pfe1, pfe2, pfe3, pfe4 

    PF0 = [0 for _ in range(noc)] 
    PF04cor = [0 for _ in range(noc)] 

    #from scipy.special import binom
    #binom(4,1) = 4 
    #binom(4,2) = 6 
    #binom(4,3) = 4 
    #binom(4,4) = 1 

    for i in range(noc-1): 
    
        PF0[i] = 1 - np.power(1 - pfe0[i], Nc) 

        PF04cor[i] = 16*pfe0[i] - (4*pfe1[i] + 6*pfe2[i] + 4*pfe3[i] + pfe4[i]) - (6*pfe2[i] + 4*pfe3[i] + pfe4[i]) - (4*pfe3[i] + pfe4[i]) - pfe4[i] 
     
        PF04cor[i] = 1 - np.power(1 - PF04cor[i], Nc) 

        #print PF0[i], PF04cor[i]


    print(noc, PF04cor[0])


if __name__ == "__main__":
    main()
