import numpy as np
from scipy.optimize import minimize
exec(open("check_core.py").read().split("# ---------- 1.")[0])   # réutilise les briques
rng = np.random.default_rng(7)

def whittle_obj(par, x, p, q):
    d = par[0]; ar = par[1:1+p]; ma = par[1+p:1+p+q]
    n = len(x); m = (n-1)//2; j = np.arange(1, m+1); lam = 2*np.pi*j/n
    I = np.abs(np.fft.fft(x)[1:m+1])**2/n
    z = np.exp(-1j*lam)
    th = 1 + sum(ma[i]*z**(i+1) for i in range(q)); ph = 1 - sum(ar[i]*z**(i+1) for i in range(p))
    g = np.abs(th)**2/np.abs(ph)**2*(2*np.sin(lam/2))**(-2*d)
    return np.log(np.mean(I/g)) + np.mean(np.log(g))


def dl_full(x, gam):    # version avec vecteur de uns (profilage GLS de la moyenne)
    n=len(x); v=np.empty(n); e=np.empty(n); e1=np.empty(n)
    v[0]=gam[0]; e[0]=x[0]; e1[0]=1.0; phi=np.array([])
    for k in range(1,n):
        kap = gam[1]/v[0] if k==1 else (gam[k]-np.sum(phi*gam[k-1:0:-1]))/v[k-1]
        phi=np.r_[phi-kap*phi[::-1],kap]; v[k]=v[k-1]*(1-kap**2)
        e[k]=x[k]-np.sum(phi*x[k-1::-1][:k]); e1[k]=1-np.sum(phi)
    return e,e1,v

# 1. moyenne GLS via innovations == GLS direct
n=120; gam=farima_acov(n,0.3,np.array([0.3]),np.array([]))
from scipy.linalg import toeplitz
G=toeplitz(gam); x=np.linalg.cholesky(G)@rng.standard_normal(n)+0.7
e,e1,v=dl_full(x,gam); mu=np.sum(e1*e/v)/np.sum(e1**2/v)
Gi=np.linalg.inv(G); one=np.ones(n); mu_d=(one@Gi@x)/(one@Gi@one)
print("mu GLS innovations = %.8f ; direct = %.8f" % (mu, mu_d))

# 2. Davies-Harte : autocov empirique vs théorique
def davies_harte(n,d,ar,ma):
    gam=farima_acov(n+1,d,ar,ma); m=2*n
    c=np.r_[gam[:n+1],gam[n-1:0:-1]]; lam=np.real(np.fft.fft(c))
    assert lam.min()>-1e-8, lam.min()
    w=rng.standard_normal(m)+1j*rng.standard_normal(m)
    y=np.fft.fft(np.sqrt(np.maximum(lam,0))*w)/np.sqrt(m)
    return y.real[:n]
n=512; ar,ma=np.array([0.4]),np.array([-0.2]); gt=farima_acov(n,0.3,ar,ma)
S=np.array([davies_harte(n,0.3,ar,ma) for _ in range(3000)])
emp=[np.mean(S[:,:n-h]*S[:,h:]) for h in (0,1,5,50)]
print("DH autocov emp:",np.round(emp,3)," théo:",np.round(gt[[0,1,5,50]],3))

# 3. SE asymptotique de Whittle (D) vs sd empirique, FARIMA(1,d,0)
def u2par(u,p,q):
    return u
def jac(eta,lam,p,q,h=1e-5):
    def logg(e):
        z=np.exp(-1j*lam); ph=1-sum(e[1+i]*z**(i+1) for i in range(p)); th=1+sum(e[1+p+i]*z**(i+1) for i in range(q))
        return np.log(np.abs(th)**2/np.abs(ph)**2*(2*np.sin(lam/2))**(-2*e[0]))
    return np.column_stack([(logg(eta+h*np.eye(len(eta))[i])-logg(eta-h*np.eye(len(eta))[i]))/(2*h) for i in range(len(eta))])
n=2500; m=(n-1)//2; lam=2*np.pi*np.arange(1,m+1)/n
eta0=np.array([0.2,0.5]); J=jac(eta0,lam,1,0); D=J.T@J/m; se_th=np.sqrt(np.diag(2*np.linalg.inv(D)/n))
est=[]
for r in range(60):
    x=davies_harte(n,0.2,np.array([0.5]),np.array([]))
    res=minimize(whittle_obj,[0.1,0.1],args=(x,1,0),method="L-BFGS-B",bounds=[(-.45,.45),(-.95,.95)]); est.append(res.x)
est=np.array(est); print("SE théorique:",se_th.round(4)," sd empirique:",est.std(0).round(4)," moy:",est.mean(0).round(3))
