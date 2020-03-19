def solve(P0,rStop,N):
    dr = rStop/float(N-1)
    beta=bet(E[0])
    #beta=bet(E) #arbitrary case
    eden=eos(P0)
    yh = np.zeros([N,2])  #yh = newtonian
    yh[0,0] = P0
    yh[0,1] = beta*dr**3*eden/3.
    yt = np.zeros([N,2])  #yt = tov/GR effect
    yt[0,0] = P0 
    yt[0,1] = beta*dr**3*eden/3.
    r=np.linspace(0.1,rStop,N)
    for j in range(N-1):
        yh[j+1] = rungkutta4(yh[j],r[j],dr,hydros)
        yt[j+1] = rungkutta4(yt[j],r[j],dr,tov)
    return r, yh[:,0], yh[:,1], yt[:,0], yt[:,1]

def mr_solve(P0,rStop,N):
    dr = rStop/float(N-1)
    beta=bet(E)
    eden=eos(P0)
    yh = np.zeros([N,2])
    yh[0,0] = P0
    yh[0,1] = beta*dr**3*eden/3.
    yt = np.zeros([N,2])
    yt[0,0] = P0
    yt[0,1] = beta*dr**3*eden/3.
    r=np.linspace(0.1,rStop,N)
    for j in range(N-1):
        yh[j+1] = rungkutta4(yh[j],r[j],dr,hydros)
        yt[j+1] = rungkutta4(yt[j],r[j],dr,tov)
    rstar = 0. 
    mstar = 0.
    count = 0
    for i in yh[:,0]:
        if i>1.e-20:
            count = count + 1    
            rstar = rstar + dr
            mstar = yh[count,1]
    rstar2 = 0. 
    mstar2 = 0.
    count2 = 0
    for i in yt[:,0]:
        if i>1.e-20:
            count2 = count2 + 1    
            rstar2 = rstar2 + dr
            mstar2 = yt[count2,1]
    return rstar,mstar,rstar2,mstar2

def mass_radius(pmin,pmax): 
    imax = 30    
    pc = np.zeros(imax)
    mass = np.zeros(imax)
    radius = np.zeros(imax)
    mass2 = np.zeros(imax)
    radius2 = np.zeros(imax)
    for i in range(imax):
        pc[i] = pmin + (pmax-pmin)*i/30
        radius[i],mass[i],radius2[i],mass2[i] = mr_solve(pc[i],40,500)
    return pc,radius,mass,radius2,mass2
