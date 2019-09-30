import numpy as np
import pandas as pd
import random as rn
import math as ma
import matplotlib as mpl
import matplotlib.pyplot as plt
from random import random
from scipy.spatial import distance_matrix as dm
from copy import deepcopy as dc


class peptide:
    def __init__(self, peptide=""):
        self.pst=peptide
        self.opm={
            (1,0):"E",
            (0,1):"N",
            (-1,0):"W",
            (0,-1):"S"
        }
        self.mpo=dict([v,k] for k,v in self.opm.items())

        # Initialize the peptide, by aligning residues along the X axis
        self.pep = self.pep_init(self.pst)

        # Initialize orientations
        self.ori = self.calc_ori(self.pep)

        # Initialize distance matrix
        m = [self.pep[i][1] for i in range(0,len(self.pep))]
        self.dst = dm(m,m)

    def pep_init(self, peptide):
        print(peptide)
        pep=[]
        for r in range(0,len(peptide)):
            pep.append([peptide[r],[r, 0]])
        return(pep)

    def valid(self, pep, ori):
        for i in range(0,len(pep)):
            for j in range(i+1, len(pep)):
                # print("A:", self.pep[i][1])
                # print("B:", self.pep[j][1])
                if pep[i][1] == pep[j][1]:
                    return(False)

        for k in range(0, len(ori)):
            x = ori[k]
            if x[0] != 0 and x[1] != 0:
              return(False)

        return(True)

    def calc_ori(self, pep):
        """
        Given a list of points in space return the list of directional
        changes from one point to the next.
        """
        ori=[]
        for r in range(1, len(pep)):
           p = pep[r-1][1]
           q = pep[r][1]
           ori.append([q[i] - p[i] for i in range(0,len(p))])

        return(ori)

    def get_ori_string(self):
        ors = []
        for i in self.ori():
            ors.append(self.mpo[tuple(self.ori[i])])

    def op_update(self, pep, ori):
        """ 
        Verify that the proposed configuration is geometrically valid.
        If the configuration is valid, then update the status of the
        peptide.
        """

        if self.valid(pep, ori):
            self.pep = pep
            self.ori = ori
            m = [self.pep[i][1] for i in range(0,len(self.pep))]
            self.dst = dm(m,m)
            return(True)
        else:
            return(False)

    def write_peptide(self, moves):
        """
        the "moves" argument expects a list of directional movements, one less
        the length of the peptide string.  All movements reference
        a global X, Y coordinate system: values in X increase to the east and
        decrease to the west.  Values in Y increase to the north and decrease
        to the south. The first residue of the peptide is anchored at position
        0,0.  All subsequent reside positions are defined by integer
        coordinates.
        """
        pep = dc(self.pep)

        if len(pep) - len(moves) != 1:
            print("Invalid move list")
            return("ER3")

        for i in range(1,len(moves)):
            po = moves[i-1]
            m = moves[i]
            if ((po == 'N' and m == 'S') or 
                (po == 'E' and m == 'W') or
                (po == 'S' and m == 'N') or
                (po == 'W' and m == 'E')):
                print("Immediate Backtracking is Invalid")
                return("ER4")

        cp = [0,0]
        pep[0][1] = cp
        for j in range(1, len(pep)):
            if moves[j-1] == "N":
                cp = [cp[i] + [0,1][i] for i in range(0,len(cp))]
                pep[j][1] = cp

            if moves[j-1] == "S":
                cp = [cp[i] + [0,-1][i] for i in range(0,len(cp))]
                pep[j][1] = cp

            if moves[j-1] == "E":
                cp = [cp[i] + [1,0][i] for i in range(0,len(cp))]
                pep[j][1] = cp

            if moves[j-1] == "W":
                cp = [cp[i] + [1,0][i] for i in range(0,len(cp))]
                pep[j][1] = cp
        return(self.op_update(pep, self.calc_ori(pep)))
    
    def write_config(self):
        c = []
        for i in self.ori:
            c.append(self.opm[tuple(i)])
        sep=""
        return(sep.join(c))
 
    def clock(self, r):
        """
        perform a clockwise rotation of the residue about its upstream 
        neighbor; apply drag transition to downstream residues
        """
        pep=dc(self.pep)
        if r < 1:
            print("Cannont Move Residude 0")
            return "ER1"
        if self.ori[r-1] == [1, 0]:
            #print("South")
            for i in range(r, len(self.pep)):
              cp = pep[i][1]
              mv = [-1, -1]
              pep[i][1] =  [cp[j] + mv[j] for j in range(0, len(cp))]
              ori = self.calc_ori(pep)

        elif self.ori[r-1] == [0, -1]:
            #print("West")
            for i in range(r, len(self.pep)):
              cp = pep[i][1]
              mv = [-1, 1]
              pep[i][1] =  [cp[j] + mv[j] for j in range(0, len(cp))]
              ori = self.calc_ori(pep)

        elif self.ori[r-1] == [-1, 0]:
            #print("North")
            for i in range(r, len(self.pep)):
              cp = pep[i][1]
              mv = [1, 1]
              pep[i][1] =  [cp[j] + mv[j] for j in range(0, len(cp))]
              ori = self.calc_ori(pep)

        elif self.ori[r-1] == [0, 1]:
            #print("East")
            for i in range(r, len(self.pep)):
              cp = pep[i][1]
              mv = [1, -1]
              pep[i][1] =  [cp[j] + mv[j] for j in range(0, len(cp))]
              ori = self.calc_ori(pep)

        else:
            return("ER2")
        return(self.op_update(pep, ori))

    def counter(self, r):
        """
        perform a counter clockwise rotation of the residu
        about its upstream neighbor; apply drag transition
        to downstream residues.
        """
        pep=dc(self.pep)
        if r < 1:
          return "ER1"
        if self.ori[r-1] == [1, 0]:
            #print("North")
            for i in range(r, len(self.pep)):
              cp = pep[i][1]
              mv = [-1, 1]
              pep[i][1] =  [cp[j] + mv[j] for j in range(0, len(cp))]
              ori = self.calc_ori(pep)

        elif self.ori[r-1] == [0, 1]:
            #print("West")
            for i in range(r, len(self.pep)):
              cp = self.pep[i][1]
              mv = [-1, -1]
              pep[i][1] =  [cp[j] + mv[j] for j in range(0, len(cp))]
              ori = self.calc_ori(pep)

        elif self.ori[r-1] == [-1, 0]:
            #print("South")
            for i in range(r, len(self.pep)):
              cp = self.pep[i][1]
              mv = [1, -1]
              pep[i][1] =  [cp[j] + mv[j] for j in range(0, len(cp))]
              ori = self.calc_ori(pep)

        elif self.ori[r-1] == [0, -1]:
            #print("East")
            for i in range(r, len(self.pep)):
              cp = self.pep[i][1]
              mv = [1, 1]
              pep[i][1] =  [cp[j] + mv[j] for j in range(0, len(cp))]
              ori = self.calc_ori(pep)

        else:
            return("ER2")
        return(self.op_update(pep, ori))

    def calc_energy(self):
        hr = []
        for i in range(0, len(self.pst)):
            if self.pst[i] == "H":
                hr.append(i)
        # print(hr)
        hdm = self.dst[hr, :]
        hdm = hdm[:, hr]
        E = 0 - (sum(hdm[np.where(hdm == 1)]) / 2)
        return(E)


class peptide_sim():
    def __init__(self, pst, maxiter=10, T=10, fname="peptide_data"):
        self.peptide = peptide(pst)
        self.maxiter = maxiter
        self.T = T
        self.E = self.peptide.calc_energy()
        self.C = self.peptide.write_config()
        self.fname = fname
        self.model = {self.C:self.E}

    def boltzman(self, en):
        pr = (ma.exp(1)**-(en / self.T)) / (ma.exp(1)**-(self.E / self.T))
        return(pr)

    def simulate(self):
        i = 0
        while i < self.maxiter:
            r = rn.randint(1, len(self.peptide.ori))
            d = rn.randint(0,1)
            if d == 1:
                v=self.peptide.clock(r)
            else:
                v=self.peptide.counter(r)
            if v:
                i = i + 1
                en = self.peptide.calc_energy()
                if en < self.E:
                    self.E = self.peptide.calc_energy()
                    self.C = self.peptide.write_config()
                    self.model[self.C]=self.E

                elif rn.random() < self.boltzman(en):
                    self.E = self.peptide.calc_energy()
                    self.C = self.peptide.write_config()
                    self.model[self.C]=self.E

    def write_model(self):
        sep="."
        fn = sep.join([self.fname, "txt"])
        data = {}
        data['Config'] = list(self.model.keys())
        data['Energy'] = list(self.model.values())

        df=pd.DataFrame.from_dict(data)
        df.to_csv(fn,index=False)

    def plot_hist(self):
        sep="."
        fn = sep.join([self.fname, "png"])
        x=list(self.model.values())
        x=[0-x[i] for i in range(0,len(x))]
        n, bins, patches = plt.hist(x, 14, density=True, facecolor='g', alpha=0.75)
        plt.xlabel('Energy ')
        plt.ylabel('Probability')
        plt.title('Distribution of Energy Values')
        plt.xlim(6, 20 )
        plt.ylim(0, 2)
        plt.grid(True)
        plt.savefig(fn)
