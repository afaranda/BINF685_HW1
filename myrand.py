# Implement my own version of random.choices based on random.choice
import random as rn

def choices(population, k=1):
    out = []
    if len(population) > 0:
        for i in range(0, k):
            out.append(rn.choice(population))
    return(out)
