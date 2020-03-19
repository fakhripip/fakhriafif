from scipy.constants import G,c
Msun=1.989e30 #kilogram
R0=G*Msun/(1.0e3*c**2)
    
def convert(eps0):
    return eps0*1.0e9/1.78762405050753e+47

def eps(a,K,gamma):
    f=((R0/a)**gamma/K)**(1./(gamma-1))
    f1=convert(f/10)
    return f1,f

def kbar(K,gamma):
    e0=E[1]  #e0
    return K*e0**(gamma-1)

def bet(e0):
    return 4*pi*e0


#EOSnya:
def eos_nr(P):
    K=kbar(knr,gnr)   #kbar
    if P<=0:
        P=1.0e-20
    return (P/K)**(3./5.)

def eos_r(P):
    K=1./3.
    if P<=0:
        P=1.0e-20
    return (P/K)**(3./4.)

def eos(P):
    # this is the fit to the non-interacting neutron matter pressure 
    anr = 2.53784142
    ar = 2.79036279
    # this is the fit to the non-interacting neutron, proton, & electron matter pressure 
    #anr = 2.33495528
    #ar = 2.94745935
    # for nuclear interaction
    #A = 0.87244601
    if P<=0: 
        P=1.e-20
    return anr*P**(3./5.) + ar*P

#E=eps(1,knr,gnr) #limiting case
#E=convert(ep0) #arbitrary case
#----------------------------------------------------------

