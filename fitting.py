#------------------------non-interacting gas fermi--------------------------
from sympy import *

u=Symbol('u')

#electron
e_n=integrate(sqrt(u**2+1.0)*u**2) #integrasi pers. (10)
e_n=integrate(u**4/sqrt(u**2+1.0)) #integrasi pers. (13)
y1=simplify(e_n)     #simplifikasi hasil integral
y2=simplify(e_n)  #simplifikasi dengan menambah faktor 1/3 di depan persamaan (13)

#proton
e_p=integrate(sqrt(u**2+(m_p/m_n)**2)*u**2)
p_p=integrate(u**4/sqrt(u**2+(m_p/m_n)**2))
y3=simplify(e_p)
y4=simplify(p_p)

#elektron
e_e=integrate(sqrt(u**2+(m_e/m_n)**2)*u**2)
p_e=integrate(u**4/sqrt(u**2+(m_e/m_n)**2))
y5=simplify(e_e)
y6=simplify(p_e)

print(y1)
print(y2)
print(y3)
print(y4)
print(y5)
print(y6)

#------------------------------

u=np.linspace(0.0,2.,100)

#neutron
ep_n=3*(0.25*u**5 + 0.375*u**3 + 0.125*u - 0.125*np.sqrt(1.0*u**2 + 1)*np.arcsinh(1.0*u))/np.sqrt(1.0*u**2 + 1)
pr_n=(0.25*u**5 - 0.125*u**3 - 0.375*u + 0.375*np.sqrt(1.0*u**2 + 1)*np.arcsinh(1.0*u))/np.sqrt(1.0*u**2 + 1)
#proton
ep_p=3*(0.250344604629826*u**5 + 0.374483804588575*u**3 + 0.124484514816418*u - 0.124313159255507*np.sqrt(1.00275873707623*u**2 + 1)*np.arcsinh(1.00137841851931*u))/np.sqrt(1.00275873707623*u**2 + 1)
pr_p=(0.250344604629827*u**5 - 0.124827934862858*u**3 - 0.373453544449255*u + np.sqrt(1.00275873707623*u**2 + 1)*(-4.44089209850063e-17*u**4 + 2.22044604925031e-17*u**2 + 0.372939477766521)*np.arcsinh(1.00137841851931*u))/np.sqrt(1.00275873707623*u**2 + 1)
#electron
ep_e=3*(459.670915152463*u**5 + 0.000203950254213724*u**3 + 2.01089300205833e-11*u + (1.09464914351283e-23*u**2 - 1.09365903724372e-14)*np.sqrt(3380757.60379364*u**2 + 1)*np.arcsinh(1838.68366060985*u))/np.sqrt(3380757.60379364*u**2 + 1)
pr_e=(459.670915152463*u**5 - 6.79834180712413e-5*u**3 - 6.032679006175e-11*u + (3.28097711173121e-14 - 2.11758236813575e-23*u**2)*np.sqrt(3380757.60379364*u**2 + 1)*np.arcsinh(1838.68366060985*u))/np.sqrt(3380757.60379364*u**2 + 1)

ep_total=ep_n+ep_p+ep_e
pr_total=pr_n+pr_p+pr_e

#------------------------------
plt.plot(ep_total,pr_total)
np.savetxt("p vs e_multi.txt",list(zip(ep_total,pr_total)), fmt="%12.5e")

#--------------------------

from scipy.optimize import curve_fit

def eos_N(P,ANR,AR):
    return ANR*P**(3./5.)+AR*P

e,P=np.loadtxt("p vs e_multi.txt",unpack=True)

popt,pcov=curve_fit(eos_N,P,e)

#-----------------Nuclear interaction-----------------------
EF0 = 22.1 #Mev
A = -122.1
B = 65.39
sig = 2.112
u = np.arange(0,10,0.01)
n0 = 0.16

#symmetry
BE0 = EF0*u**(2/3) + A/2*u + B/(sig+1)*u**sig
P0 = n0*(2/3*EF0*u**(5/3) + A/2*u**2 + sig*B/(sig+1)*u**(sig+1))

m_N = 939
S0 = 30

#asymmetry
def enden(alpha,k):      
    E = EF0*u**(2/3) + A/2*u + B/(sig+1)*u**sig + alpha**2*((2**(2/3)-1)*EF0*(u**(2/3) - u**k) + S0*u**k)
    return u,E

def pres(alpha,k):
    P = n0*(2/3*EF0*u**(5/3) + A/2*u**2 + sig*B/(sig+1)*u**(sig+1)) + n0*alpha**2*((2**(2/3)-1)*EF0*(2/3*u**(5/3) - k*u**(k+1)) + S0*k*u**(k+1))
    return u,P

def en(alpha,k):      
    E = m_N + EF0*u**(2/3) + A/2*u + B/(sig+1)*u**sig + alpha**2*((2**(2/3)-1)*EF0*(u**(2/3) - u**k) + S0*u**k)
    return u,E

E1=en(1,1)
P1_1=pres(1,1)
plt.plot(E1[0],E1[1])
plt.plot(P1_1[0],P1_1[1],color="red")
plt.xlim(0,10)
plt.ylim(0,10000)

#asymmetry dimensionless for fitting
ep0=m_n**4*c**5/(3.*pi**2*hbar**3)*6.25e-33 #MeV
def enden2(alpha,k):      
    E = (EF0*u**(2/3) + A/2*u + B/(sig+1)*u**sig + alpha**2*((2**(2/3)-1)*EF0*(u**(2/3) - u**k) + S0*u**k))/ep0
    return u,E

def pres2(alpha,k):
    P = (n0*(2/3*EF0*u**(5/3) + A/2*u**2 + sig*B/(sig+1)*u**(sig+1)) + 
         n0*alpha**2*((2**(2/3)-1)*EF0*(2/3*u**(5/3) - k*u**(k+1)) + S0*k*u**(k+1)))/ep0    
    return u,P

def en2(alpha,k):      
    E = (m_N + EF0*u**(2/3) + A/2*u + B/(sig+1)*u**sig + alpha**2*((2**(2/3)-1)*EF0*(u**(2/3) - u**k) + S0*u**k))/ep0
    return u,E

E2=en2(1,1)
P2_1=pres2(1,1)
#------------------------------------------------
plt.plot(P2_1[1],E2[1])
np.savetxt("P2 vs BE2.txt",list(zip(P2_1[1],E2[1])), fmt="%12.5e")
#-----------------------------------------------
from scipy.optimize import curve_fit
import numpy as np

def eos_nuc(p,k0):
    return k0*p**(1/2)

p2,e2=np.loadtxt("P2 vs BE2.txt",unpack=True)

popt,pcov=curve_fit(eos_nuc,p2,e2)
print(popt)

#-------------------------------causality----------------------------

EF0 = 22.1 #Mev
A = -122.934
B = 70.1
C = 0.15
sig = 2.
m_N = 939.
S0 = 30.
n0 = 0.16
u = np.arange(0,10,0.01)   
ep0=m_n**4*c**5/(3.*pi**2*hbar**3)*6.25e-33 #MeV 

def enden34(alpha,k):      
    E = (EF0*u**(2/3) + A/2*u + B/(sig+1)*u**sig/(1+C*u**(sig-1)) + alpha**2*((2**(2/3)-1)*EF0*(u**(2/3) - u**k) + S0*u**k))/ep0
    return u,E

def pres4(alpha,k):
    P = (n0*(2/3*EF0*u**(5/3) + A/2*u**2 + B*(C*u**(2*sig)+sig*u**(sig+1))/((sig+1)*(C*u**(sig-1)+1)**2)) + n0*alpha**2*((2**(2/3)-1)*EF0*(2/3*u**(5/3) - k*u**(k+1)) + S0*k*u**(k+1)))/ep0
    return u,P
            
def en4(alpha,k):      
    Enn = (m_N + EF0*u**(2/3) + A/2*u + B/(sig+1)*u**sig/(1+C*u**(sig-1)) + alpha**2*((2**(2/3)-1)*EF0*(u**(2/3) - u**k) + S0*u**k))/ep0
    return u,Enn

End3 =en4(1,1)
P4_1=pres4(1,1)
#-----------------------------------------------
plt.plot(P4_1[1],End3[1])
np.savetxt("P4 vs BE4.txt",list(zip(P4_1[1],End3[1])), fmt="%12.5e")
#--------------------------------------
from scipy.optimize import curve_fit
import numpy as np

def eos_nuc(p,k0):
    return k0*p**(1/2)

p4,e4=np.loadtxt("P4 vs BE4.txt",unpack=True)

popt,pcov=curve_fit(eos_nuc,p4,e4)
print(popt)

#------------finding parameters A,B, and sigma-----------------
EF0 = 22.1
BE = -16

def B(k0):
    sig = (k0 + 2*EF0) / (3*EF0 - 9*BE)
    B = (sig+1)/(sig-1)*(1/3*EF0 - BE)
    A = BE - (5/3*EF0) - B
    return A,B,sig

def y(b,sig):
    return b*(sig-1)/(sig+1)

y2=B(400)
y(70.1,2)



