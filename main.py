from darwin import darwin as dw
from pep_lat import peptide_sim as ps

s=""
pst="PPPPPPPHPPHPPPPPPPHPHPHPHHHHPHPPPPHPHPHPHPP"
pst= s.join([pst, "HPPPHPHPHPPHPPPPHPPPHPHHPPHPPHPHPHHHP"])

p=ps(pst, maxiter=1000, T=10, fname ="peptide_data_T10")
p.random_simulate()
p.write_model()

p=ps(pst, maxiter=1000, T=5, fname = "peptide_data_T5")
p.random_simulate()
p.write_model()


p=ps(pst, maxiter=1000, T=5, fname = "peptide_data_T1")
p.random_simulate()
p.write_model()

p=ps(pst, maxiter=1000, T=5, fname = "peptide_data_RN")
p.random_pop()
p.write_model()


# Initialize Genetic Algorithm
d = dw(pst, maxgen=25)

# Run the Genetic Algorithm, Print best configuration and score
d.evolution()
print(d.bst)
print(d.enr)
