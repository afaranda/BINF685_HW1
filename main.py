from pep_lat import peptide
from pep_lat import peptide_sim as ps

import random as rn


s=""
pst="PPPPPPPHPPHPPPPPPPHPHPHPHHHHPHPPPPHPHPHPHPP"
pst= s.join([pst, "HPPPHPHPHPPHPPPPHPPPHPHHPPHPPHPHPHHHP"])

p=ps(pst, maxiter=100000)
p.simulate()
p.write_model()
