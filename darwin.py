import random as rn
import math as ma
from random import random
from pep_lat import peptide_sim as ps
from pep_lat import peptide as pp
from copy import deepcopy as dc

class darwin:
    def __init__(self, pst="HPPHP", init_T = 10, fwd_carry = 100,
                popsize = 1000, ncross=1, nmut=20, maxgen = 10):
        self.pst = pst
        self.itm = init_T
        self.fwd = fwd_carry
        self.psz = popsize
        self.ncr = ncross
        self.nmt = nmut
        self.pop = None
        self.bst = None
        self.enr = None
        self.gen = maxgen
        self.peptide = pp(pst)
    
    def make_pop(self):
        p = ps(pst=self.pst, maxiter=self.psz, num=self.nmt)
        p.random_pop()
        self.pop = p.model
        self.bst = sorted(self.pop, key=self.pop.get)[0]

    def boltzman(self, en, ave, t):
        pr = (ma.exp(1)**-(en / t)) / (ma.exp(1)**-(ave / t))
        return(pr)

    def cross(self, cf, cm):
        ncr = self.ncr
        if ncr > (len(cf) // 10):
            print("Too many crossovers for this peptide")
            return(False)
        crs = [0]
        crs = crs + (rn.sample(range(1,len(cf) - 1), ncr))
        crs.append(len(cf))
        # print(crs)
        crs.sort()
        child1 = []
        child2 = []
        last = 0
        for i in range(0, ncr+2):
            if i % 2 == 1:
                child1.append(cm[last:crs[i]])
                child2.append(cf[last:crs[i]])
                last = crs[i]
            if i % 2 == 0:
                child1.append(cf[last:crs[i]])
                child2.append(cm[last:crs[i]])
                last = crs[i]
        sep=""
        child1 = sep.join(child1)
        child2 = sep.join(child2)
        return(child1, child2)

    def mutate_pop(self, pop):
        p = ps(self.pst, num = min(1, self.nmt // 3))
        newpop = {}
        for i in pop.keys():
            cf = p.mutate(i)
            p.peptide.write_peptide(cf)
            newpop[p.peptide.write_config()] = p.peptide.calc_energy()
        return(newpop)

    def run_crosses(self, pop):
        newpop={}
        pop=sorted(pop, key=pop.get)
        cf = pop.pop(0)
        n = 3 # Necessary for decay function
        t = self.itm
        while len(pop) > 0:
            cm = pop.pop(0)
            child1, child2 = self.cross(cf, cm)
            # print(child1)
            # print(child2)
            p = pp(self.pst)
            q = pp(self.pst)
            if p.write_peptide(child1) and q.write_peptide(child2):
                pe, qe = p.calc_energy(), q.calc_energy()
                ave = sum([pe, qe])/2
                for i in [p, q]:
                    ie = i.calc_energy()
                    #print(ie, ave)
                    if ie <= ave:
                        newpop[i.write_config()] = ie
                    elif rn.random() < self.boltzman(ie, ave, t):
                        # print(self.boltzman(ie, ave, t))
                        newpop[i.write_config()] = ie
                
                t = t / (ma.log(n))
                n = n + 1
                if len(pop) > 1:
                    cf = pop.pop(0) 
        return(newpop)

    def update_pop(self, pop1, pop2):
        fwd = self.fwd
        opk = sorted(pop1, key=pop1.get)[0:(len(pop1) - fwd)]
        op = {key:pop1[key] for key in opk}
        npk = set(pop2.keys()).difference(set(pop1.keys()))
        np = {key:pop2[key] for key in npk}
        npk = sorted(np, key=np.get)[0:fwd]
        np = {key:np[key] for key in npk}
        op.update(np)
        return(op)  

    def evolution(self):
        self.make_pop()
        # lbs = self.best
        g = 0
        while g < self.gen:
            pop = dc(self.pop)
            pop2 = self.mutate_pop(pop)
            pop2 = self.run_crosses(pop2)
            self.pop=self.update_pop(pop, pop2)

            p=pp(self.pst)
            p.write_peptide(self.bst)
            
            q=pp(self.pst)
            q.write_peptide(sorted(self.pop, key=self.pop.get)[0])
            if q.calc_energy() < p.calc_energy():
                self.bst = q.write_config()
                self.enr = q.calc_energy()
            g = g + 1    
                
        


