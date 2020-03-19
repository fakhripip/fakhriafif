#ans=solve(1.0e-4,40,500)

#untuk print nilai M & R
print("saat input P(0) = 0.05")
print("beta = {}".format(bet(E))
print("eps0 = {}".format(E))
print("------Newtonian-------")
print("R = {}".format(ans[0]))
print("M = {}".format(ans[1]))
print("-------GR-------------")
print("R = {}".format(ans[2]))
print("M = {}".format(ans[3]))

untuk plotting P-R & M-R

fig=plt.figure(figsize=(9.5,6.5))
ax1=fig.add_subplot(211)
ax1.plot(ans[0],ans[1],color='blue',label='Newtonian')
ax1.plot(ans[0],ans[3],color='red',label='GR') #GR
ax1.set_xlabel('radius(km)')
ax1.set_ylabel(r'Pressure(ergs/cm$^3$)')
ax1.legend(loc='upper right')
ax1.set_xlim(0,17)
ax1.set_ylim(0,1.0e-4)

ax2=fig.add_subplot(212)
ax2.plot(ans[0],ans[2],color='blue',label='Newtonian')
ax2.plot(ans[0],ans[4],color='red',label='GR')
ax2.set_xlabel('radius(km)')
ax2.set_ylabel(r'Massa (M$_\odot$)')
ax2.legend(loc='lower right')
ax2.set_xlim(0,17)
ax2.set_ylim(0,0.8)

fig.tight_layout()
#fig.savefig('nonrel(1e-4).png')
fig.show()

#plotting relasi M_R sebagai fungsi P(0)
#ans=mass_radius(0.001,5.5) 

fig=plt.figure(figsize=(6,4))
plt.plot(ans[3],ans[4],color='red',label='GR') #GR
plt.xlabel('radius(km)')
plt.ylabel(r'Massa (M$_\odot$)')
plt.legend(loc='upper right')
#plt.xlim(4,15)
#plt.ylim(0.5,0.9)

fig.tight_layout()
#fig.savefig('comb GR(M-R relation).png')
fig.show()
