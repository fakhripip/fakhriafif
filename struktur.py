from scipy.constants import G,c
Msun=1.989e30 #kilogram
R0=G*Msun/(1.0e3*c**2)

def hydros(state,r):
    beta=bet(E[0])
    #beta=bet(E) #arbitrary case
    eden=eos(state[0])
    one=-R0*eden*state[1]/r**2
    two=beta*r**2*eden
    return np.array([one,two])

def tov(state,r):
    beta=bet(E[0])
    #beta=bet(E) #arbitrary case
    eden=eos(state[0])
    one=-R0*eden*state[1]*(1+(state[0]/eden))*(1 + beta*(state[0]/state[1])*r**3)/(r**2-2.*R0*state[1]*r)
    two=beta*r**2*eden
    return np.array([one,two])
