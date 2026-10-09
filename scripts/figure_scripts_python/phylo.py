import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family']='Liberation Sans'; plt.rcParams['mathtext.fontset']='custom'; plt.rcParams['mathtext.it']='Liberation Sans:italic'; plt.rcParams['mathtext.bf']='Liberation Sans:bold:italic'
fig,ax=plt.subplots(figsize=(6.3,3.5))
# tips: y, label, symbol, study
tips={1:(r'$\mathbf{Exaiptasia\ diaphana}$ (sea anemone)','dino',True),
      2:('Reef-building corals','dino',False),
      3:(r'$\it{Hydra\ viridissima}$ (hydrozoan)','green',False),
      4:(r'$\mathbf{Waminoa}$ sp. (acoel flatworm)','dino',True),
      5:(r'Giant clams ($\it{Tridacna}$)','dino',False),
      6:('Deuterostomes (e.g. sea urchins, vertebrates)',None,False)}
y={k:7-k for k in tips}  # top to bottom
xt=5.0
L=dict(color='black',lw=1.2)
def h(x0,x1,yy,**k): ax.plot([x0,x1],[yy,yy],**{**L,**k})
def v(x,y0,y1,**k): ax.plot([x,x],[y0,y1],**{**L,**k})
# Anthozoa node
xa=3.6; ya=(y[1]+y[2])/2; h(xa,xt,y[1]); h(xa,xt,y[2]); v(xa,y[2],y[1])
# Cnidaria node
xc=2.2; yc=(ya+y[3])/2; h(xc,xa,ya); h(xc,xt,y[3]); v(xc,y[3],ya)
# Nephrozoa node
xn=3.6; yn=(y[5]+y[6])/2; h(xn,xt,y[5]); h(xn,xt,y[6]); v(xn,y[6],y[5])
# Bilateria node
xb=2.2; yb=(y[4]+yn)/2; h(xb,xt,y[4]); h(xb,xn,yn); v(xb,yn,y[4])
# root
xr=0.8; yr=(yc+yb)/2; h(xr,xc,yc); h(xr,xb,yb); v(xr,yb,yc); h(0.4,xr,yr)
lab=dict(fontsize=8,ha='right',va='bottom',color='0.25')
ax.text(xa-0.05,ya+0.08,'Anthozoa',**lab); ax.text(xc-0.05,yc+0.08,'Cnidaria',**lab)
ax.text(xn-0.05,yn+0.08,'Nephrozoa',**lab); ax.text(xb-0.05,yb+0.08,'Bilateria',**lab)
ax.text(xt-0.05,y[4]+0.08,'Xenacoelomorpha',**lab)
for k,(t,s,st) in tips.items():
    if s=='dino': ax.plot(xt+0.18,y[k],'o',ms=7,mfc='#C8813A',mec='black',mew=0.6)
    elif s=='green': ax.plot(xt+0.18,y[k],'^',ms=7.5,mfc='#4C9A2A',mec='black',mew=0.6)
    ax.text(xt+0.42,y[k],t,fontsize=9,va='center')
ax.plot([],[],'o',ms=7,mfc='#C8813A',mec='black',mew=0.6,ls='',label='dinoflagellate symbionts')
ax.plot([],[],'^',ms=7.5,mfc='#4C9A2A',mec='black',mew=0.6,ls='',label='green-algal symbionts')
ax.legend(loc='lower left',bbox_to_anchor=(0.0,-0.13),frameon=False,fontsize=8,ncol=2,handletextpad=0.3)
ax.set_xlim(0.3,10.3); ax.set_ylim(0.2,6.6); ax.axis('off')
fig.tight_layout(); fig.savefig('fig03_phylogeny.png',dpi=600)
