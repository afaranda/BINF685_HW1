from darwin import darwin as dw
from pep_lat import peptide_sim as ps

s=""
pst="PPPPPPPHPPHPPPPPPPHPHPHPHHHHPHPPPPHPHPHPHPP"
pst= s.join([pst, "HPPPHPHPHPPHPPPPHPPPHPHHPPHPPHPHPHHHP"])

p=ps(pst, maxiter=1000, T=10, fname ="peptide_data_T10.txt")
p.random_simulate()
p.write_model()

p=ps(pst, maxiter=1000, T=5, fname = "peptide_data_T5.txt")
p.random_simulate()
p.write_model()


p=ps(pst, maxiter=1000, T=5, fname = "peptide_data_T1.txt")
p.random_simulate()
p.write_model()



d = dw(pst)
d.evolution()
print(d.bst)
print(d.enr)
