def rungkutta2(y,time,dt,derivs):
    k0=dt*derivs(y,time)
    k1=dt*derivs(y+k0,time+dt)
    y_next=y+0.5*(k0+k1)
    return y_next

def rungkutta4(y,time,dt,derivs):
    k0=dt*derivs(y,time)
    k1=dt*derivs(y+k0/2,time+dt/2)
    k2=dt*derivs(y+k1/2,time+dt/2)
    k3=dt*derivs(y+k2,time+dt)
    y_next=y+(k0+2*k1+2*k2+k3)/6
    return y_next

def euler1(y,t,dt,derivs):
    y_next=y+derivs(y,t)*dt
    return y_next

def factorial(n):
    f=1
    for i in range(2,n+1):
        f=f*i
    return f

def sq(x):
    x=x*x
    return x
    
def root_bisection(f,a,b,tolerance=1e-6):
    dx=abs(b-a)
    while dx>tolerance:
        x=(a+b)/2.0
        if(f(a)*f(x))<0:
            b=x
        else:
            a=x
        dx=abs(b-a)
    return x

def root_newton(f,df,guess,tolerance=1e-6):
    dx=2*tolerance
    while dx>tolerance:
        x1=x-f(x)/df(x)
        dx=abs(x-x1)
        x=x1
    return x