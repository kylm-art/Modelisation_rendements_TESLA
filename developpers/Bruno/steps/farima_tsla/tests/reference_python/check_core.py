import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad
from scipy.special import gammaln
from scipy.signal import lfilter
from scipy.stats import multivariate_normal
from scipy.linalg import toeplitz, solve
rng = np.random.default_rng(1)

# ---------- briques (miroir exact de ce qui sera écrit en R) ----------
def fi_acov(d, H):                      # gamma_U(0..H), sigma^2 = 1
    g = np.empty(H+1); g[0] = np.exp(gammaln(1-2*d) - 2*gammaln(1-d))
    for h in range(1, H+1): g[h] = g[h-1]*(h-1+d)/(h-d)
    return g

def arma_psi(ar, ma, L):                # psi_0..psi_L  (ARMAtoMA + 1)
    ar, ma = np.atleast_1d(ar), np.atleast_1d(ma)
    num = np.r_[1, ma]; den = np.r_[1, -ar]
    imp = np.zeros(L+1); imp[0] = 1
    return lfilter(num, den, imp)

def arma_acov(ar, ma, M):               # gamma_ARMA(0..M), sigma^2=1 via psi (troncature)
    psi = arma_psi(ar, ma, 5000)
    return np.array([np.sum(psi[:len(psi)-h]*psi[h:]) for h in range(M+1)])

def farima_acov(n, d, ar, ma):
    ar = np.atleast_1d(ar) if len(np.atleast_1d(ar)) else np.array([]); ma = np.atleast_1d(ma) if len(np.atleast_1d(ma)) else np.array([])
    if len(ar)==0 and len(ma)==0: return fi_acov(d, n-1)
    M = 400
    a = arma_acov(ar, ma, M)
    gu = fi_acov(d, n-1+M)
    out = np.empty(n)
    hh = np.arange(-M, M+1); ah = a[np.abs(hh)]
    for k in range(n): out[k] = np.sum(ah * gu[np.abs(k-hh)])
    return out

def dl_nll(x, gam):                     # -2 loglik concentre (sigma^2 profile)
    n = len(x); v = np.empty(n); e = np.empty(n)
    v[0] = gam[0]; e[0] = x[0]; phi = np.array([])
    for k in range(1, n):
        if k == 1: kap = gam[1]/v[0]
        else:      kap = (gam[k] - np.sum(phi*gam[k-1:0:-1]))/v[k-1]
        phi = np.r_[phi - kap*phi[::-1], kap]
        v[k] = v[k-1]*(1-kap**2)
        e[k] = x[k] - np.sum(phi*x[k-1::-1][:k])
    s2 = np.mean(e**2/v)
    return n*np.log(s2) + np.sum(np.log(v)), s2, e, v

# ---------- 1. FI(d) vs intégrale spectrale ----------
d0 = 0.3
for h in (0,1,5,20):
    f = lambda l: (2*np.sin(l/2))**(-2*d0)*np.cos(h*l)/(2*np.pi)
    num = 2*quad(f, 0, np.pi, limit=400)[0]
    print("FI acov h=%2d  formule=%.8f  integrale=%.8f" % (h, fi_acov(d0,h)[h], num))

# ---------- 2. FARIMA(1,d,1) convolution vs intégrale spectrale ----------
d0, ar, ma = 0.25, np.array([0.5]), np.array([0.3])
gam = farima_acov(60, d0, ar, ma)
def spec(l):
    z = np.exp(-1j*l)
    return abs(1+ma[0]*z)**2/abs(1-ar[0]*z)**2*(2*np.sin(l/2))**(-2*d0)/(2*np.pi)
for h in (0,1,10,50):
    num = 2*quad(lambda l: spec(l)*np.cos(h*l), 0, np.pi, limit=800, points=[1e-6])[0]
    print("FARIMA acov h=%2d conv=%.7f integ=%.7f" % (h, gam[h], num))

# ---------- 3. vraisemblance Durbin-Levinson vs gaussienne directe ----------
n = 80; G = toeplitz(farima_acov(n, d0, ar, ma))
x = np.linalg.cholesky(G) @ rng.standard_normal(n)
nll, s2, e, v = dl_nll(x, farima_acov(n, d0, ar, ma))
direct = -2*multivariate_normal(mean=np.zeros(n), cov=G*s2).logpdf(x)
print("DL: -2ll = %.6f ; direct = %.6f ; s2=%.5f" % (nll + n*np.log(2*np.pi) + n, direct, s2))

# ---------- 4. simulation + récupération (Whittle et MLE exact) ----------
def sim_farima(n, d, ar, ma, sig=1.0, burn=2000, eps=None):
    N = n + burn
    e = rng.standard_normal(N)*sig if eps is None else eps
    psi_arma = arma_psi(ar, ma, N-1)
    pi_ = np.empty(N); pi_[0] = 1
    for j in range(1, N): pi_[j] = pi_[j-1]*(j-1+d)/j      # (1-B)^{-d}
    psi = np.convolve(psi_arma, pi_)[:N]
    x = np.convolve(e, psi)[:N]
    return x[burn:]

def whittle_obj(par, x, p, q):
    d = par[0]; ar = par[1:1+p]; ma = par[1+p:1+p+q]
    n = len(x); m = (n-1)//2; j = np.arange(1, m+1); lam = 2*np.pi*j/n
    I = np.abs(np.fft.fft(x - 0)[1:m+1])**2/n
    z = np.exp(-1j*lam)
    th = 1 + sum(ma[i]*z**(i+1) for i in range(q)); ph = 1 - sum(ar[i]*z**(i+1) for i in range(p))
    g = np.abs(th)**2/np.abs(ph)**2*(2*np.sin(lam/2))**(-2*d)
    return np.log(np.mean(I/g)) + np.mean(np.log(g))

n = 3000; true = dict(d=0.2, ar=np.array([0.4]), ma=np.array([-0.2]))
est_w = []
for r in range(30):
    x = sim_farima(n, true['d'], true['ar'], true['ma'])
    res = minimize(whittle_obj, [0.05,0.1,0.0], args=(x,1,1), method="L-BFGS-B",
                   bounds=[(-0.45,0.45),(-0.95,0.95),(-0.95,0.95)])
    est_w.append(res.x)
est_w = np.array(est_w)
print("Whittle  moyenne =", est_w.mean(0).round(3), " (vrai 0.2, 0.4, -0.2)  sd =", est_w.std(0).round(3))
print("SE théorique d (ARFIMA(0,d,0)) = %.4f" % np.sqrt(6/(np.pi**2*n)))

# MLE exact sur une réalisation n=600
n = 600; x = sim_farima(n, 0.2, true['ar'], true['ma'])
def exact_obj(par):
    d, a, b = par
    return dl_nll(x, farima_acov(n, d, np.array([a]), np.array([b])))[0]
r = minimize(exact_obj, [0.05,0.1,0.0], method="L-BFGS-B", bounds=[(-0.45,0.45),(-0.95,0.95),(-0.95,0.95)])
rw = minimize(whittle_obj, [0.05,0.1,0.0], args=(x,1,1), method="L-BFGS-B", bounds=[(-0.45,0.45),(-0.95,0.95),(-0.95,0.95)])
print("n=600  MLE exact =", r.x.round(3), " Whittle =", rw.x.round(3))

# ---------- 5. prévision AR(inf) tronquée vs BLP exact ----------
n = 800; d0, a0, b0 = 0.2, np.array([0.4]), np.array([-0.2])
x = sim_farima(n, d0, a0, b0); H = 10
# AR(inf): pi = Phi(B)(1-B)^d / Theta(B)
fd = np.empty(n+1); fd[0]=1
for j in range(1,n+1): fd[j] = fd[j-1]*(j-1-d0)/j
u = np.convolve(fd, np.r_[1,-a0])[:n+1]
pi_ = lfilter([1], np.r_[1,b0], u)
hist = list(x)
for h in range(H):
    K = min(len(hist), n)
    xh = np.array(hist[::-1][:K])
    hist.append(-np.sum(pi_[1:K+1]*xh))
fc_ar = np.array(hist[n:])
gam = farima_acov(n+H, d0, a0, b0)
G = toeplitz(gam[:n]); fc_ex = np.array([gam[np.arange(h+1, h+1+n)][::-1] @ solve(G, x) for h in range(H)])
print("prev. AR(inf) :", fc_ar[:5].round(4)); print("prev. exacte  :", fc_ex[:5].round(4))
print("écart max     : %.2e" % np.max(np.abs(fc_ar-fc_ex)))
# variance de prévision
psi = np.convolve(arma_psi(a0,b0,H), np.cumprod(np.r_[1,(np.arange(H)+d0)/(np.arange(1,H+1))]))[:H]
mse_h = np.cumsum(psi**2)
mse_ex = [gam[0] - gam[np.arange(h+1,h+1+n)] @ solve(G, gam[np.arange(h+1,h+1+n)]) for h in range(H)]
print("MSE psi  :", mse_h[:5].round(4)); print("MSE exact:", np.round(mse_ex[:5],4))
