#!/usr/bin/env python3
"""Dense-diagonalization plane-wave DFT oracle: Si atom in a cubic box.

Same conventions as cpu::pw: psi = (1/sqrt(Om)) sum c_G e^{iGr}, jellium
G=0 cancellations (v_H(0)=0, V_loc(0)=alpha-Z finite part, Ewald + QE G=0
term), LDA = Slater x + PW92 c, KB projectors from the UPF, integer
occupation (2 bands x 2). Gamma only.
"""
import math, re, sys
import numpy as np

UPF = 'benchmarks/pw/si/validation/reference/Si_lda_nc.upf'
L_ANG = float(sys.argv[1]) if len(sys.argv) > 1 else 8.0
ECUT_RY = float(sys.argv[2]) if len(sys.argv) > 2 else 16.0
A2B = 1.8897259886
Z = 4.0
NB = 2

txt = open(UPF).read()
def nums(tag):
    m = re.search(r'<'+tag+r'[^>]*>(.*?)</'+tag, txt, re.S)
    vals = []
    for tok in m.group(1).split():
        try: vals.append(float(tok))
        except ValueError: pass
    return vals
r = np.array(nums('PP_R')); vloc_ry = np.array(nums('PP_LOCAL'))
rhoat = np.array(nums('PP_RHOATOM'))
betas = []
for body in re.findall(r'<PP_BETA\.\d+(.*?)</PP_BETA\.\d+>', txt, re.S):
    b = []
    for tok in body.split():
        try: b.append(float(tok))
        except ValueError: pass
    betas.append(np.array(b[-len(r):]))
dij = np.array(nums('PP_DIJ')).reshape(6, 6) * 0.5   # Ry -> Ha
beta_l = [0, 0, 1, 1, 2, 2]
vloc = vloc_ry * 0.5                                  # Ry -> Ha

rm = 0.5*(r[1:]+r[:-1]); dr = np.diff(r)
vloc_m = 0.5*(vloc[1:]+vloc[:-1])
rhoat_m = 0.5*(rhoat[1:]+rhoat[:-1])
betas_m = [0.5*(bb[1:]+bb[:-1]) for bb in betas]
erf_m = np.array([math.erf(x) for x in rm])
def radint(f_mid):                                    # midpoint rule
    return float(np.sum(f_mid*dr))

L = L_ANG*A2B; om = L**3; b = 2*math.pi/L
ecut = ECUT_RY/2.0
# wfc sphere
gmax = math.sqrt(2*ecut); hmax = int(gmax/b)+1
hh = np.arange(-hmax, hmax+1)
H, K, Lm = np.meshgrid(hh, hh, hh, indexing='ij')
hkl = np.stack([H.ravel(), K.ravel(), Lm.ravel()], 1)
g2 = (hkl**2).sum(1)*b*b
sel = 0.5*g2 <= ecut+1e-12
hkl_w = hkl[sel]; g2_w = g2[sel]
order = np.argsort(g2_w, kind='stable')
hkl_w = hkl_w[order]; g2_w = g2_w[order]
ng = len(hkl_w)
print(f"wfc sphere ng={ng}")
G_w = hkl_w*b

def j0(x): return np.where(np.abs(x)<1e-8, 1-x*x/6, np.sin(np.maximum(np.abs(x),1e-30))/np.maximum(np.abs(x),1e-30))
def j1(x):
    ax = np.maximum(np.abs(x), 1e-12)
    return np.where(np.abs(x)<1e-8, x/3, np.sin(ax)/(ax*ax)-np.cos(ax)/ax)*np.sign(np.where(x==0,1,np.sign(x)))**0  # j1 odd handled via |q| use only
def j2(x):
    ax = np.maximum(np.abs(x), 1e-12)
    return np.where(np.abs(x)<1e-8, x*x/15, (3/ax**3-1/ax)*np.sin(ax)-(3/ax**2)*np.cos(ax))
jl = {0: j0, 1: j1, 2: j2}

# V_loc(G) on difference vectors via FFT grid approach:
# build V_loc(r) on grid from spherical tables? Simpler: dense matrix needs
# V(G_i - G_j): radial transform evaluated at |G_i-G_j| -> precompute by
# unique shells of the DIFFERENCE lattice (h in [-2hmax, 2hmax]).
def vloc_of_q(q):
    # q: array of |G| values -> erf-split transform (Ha), G=0 -> alpha
    out = np.empty_like(q)
    for i, G in enumerate(q):
        if G < 1e-12:
            out[i] = 4*math.pi/om * radint(rm*(rm*vloc_m+Z))
        else:
            f = rm*vloc_m + Z*erf_m
            sint = radint(rm*f*j0(G*rm)/rm)  # = int (r v + Z erf) j0 dr ... careful
            sint = radint((rm*vloc_m + Z*erf_m)*rm*j0(G*rm)/rm)
            sint = radint(rm*(rm*vloc_m)*j0(G*rm) + rm*Z*erf_m*j0(G*rm))
            out[i] = 4*math.pi/om*(sint - Z*math.exp(-G*G/4)/(G*G))
    return out

# unique |q|^2 integer shells for differences
diff_n = {}
for i in range(ng):
    d = hkl_w - hkl_w[i]
    n = (d**2).sum(1)
    for v in np.unique(n):
        diff_n[int(v)] = None
qs = np.array(sorted(diff_n.keys()))
vq = vloc_of_q(np.sqrt(qs.astype(float))*b)
vmap = dict(zip(qs.tolist(), vq.tolist()))
print("vloc shells:", len(qs), " V(0)=", vmap[0])

VL = np.zeros((ng, ng))
for i in range(ng):
    d = hkl_w - hkl_w[i]
    n = (d**2).sum(1)
    VL[:, i] = [vmap[int(x)] for x in n]   # V(G_row - G_col)

# KB projectors: columns over (proj, m)
def ylm(l, m, gv):
    x, y, z = gv[:,0], gv[:,1], gv[:,2]
    rr = np.sqrt((gv**2).sum(1))
    ir = np.where(rr>1e-16, 1/np.maximum(rr,1e-30), 0.0)
    fp = 4*math.pi
    if l==0: return np.full(len(gv), 1/math.sqrt(fp))
    if l==1:
        n1 = math.sqrt(3/fp)
        return [n1*y*ir, n1*z*ir, n1*x*ir][m]
    ir2 = ir*ir
    if m==0: return math.sqrt(15/fp)*x*y*ir2
    if m==1: return math.sqrt(15/fp)*y*z*ir2
    if m==2: return math.sqrt(5/(4*fp))*(3*z*z*ir2-1)
    if m==3: return math.sqrt(15/fp)*x*z*ir2
    return math.sqrt(15/(4*fp))*(x*x-y*y)*ir2

cols = []; Ds = []
pref = 4*math.pi/math.sqrt(om)
for ip, bet in enumerate(betas):
    l = beta_l[ip]; D = dij[ip, ip]
    bm = betas_m[ip]
    fl = np.array([radint(rm*bm*jl[l](G*rm)) for G in np.sqrt(g2_w)])
    il = {0:1, 1:-1j, 2:-1}[l]
    for m in range(2*l+1):
        cols.append(pref*il*ylm(l, m, G_w)*fl)
        Ds.append(D)
B = np.array(cols).T          # ng x ncol
Dv = np.array(Ds)
VNL = (B*Dv) @ B.conj().T

T = np.diag(0.5*g2_w)

# density: FFT grid covering 2x sphere
nfft = 2*(2*hmax+1)
nfft += nfft % 2
def to_grid(cg, hklv):
    gr = np.zeros((nfft,nfft,nfft), complex)
    idx = tuple((hklv % nfft).T)
    gr[idx] = cg
    return gr
# initial density: atomic
rho_g_at = np.zeros(ng, complex)
qsw = np.sqrt(g2_w)
for i, G in enumerate(qsw):
    rho_g_at[i] = radint(rhoat_m*j0(G*rm))/om
rho_g = rho_g_at.copy()

def pw92_exc_vxc(n):
    n = np.maximum(n, 1e-12)
    rs = (3/(4*math.pi*n))**(1/3.)
    ex = -3/(4*math.pi)*(3*math.pi**2*n)**(1/3.)
    vx = 4./3.*ex
    A,a1,b1,b2,b3,b4 = 0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294
    s = np.sqrt(rs)
    Q = 2*A*(b1*s + b2*rs + b3*rs*s + b4*rs*rs)
    ec = -2*A*(1+a1*rs)*np.log(1+1/Q)
    dQ = 2*A*(0.5*b1/s + b2 + 1.5*b3*s + 2*b4*rs)
    dec = -2*A*a1*np.log(1+1/Q) + 2*A*(1+a1*rs)*dQ/(Q*Q+Q)
    vc = ec - rs/3*dec
    return ex+ec, vx+vc

E_ewald = -2.837297*Z*Z/(2*L)    # simple cubic Madelung, point ion + background
ew_alpha_term = 0.0              # QE G=0 already inside Madelung for this form

for it in range(30):
    # V_H
    vh = np.where(g2_w>1e-12, 4*math.pi*rho_g/np.maximum(g2_w,1e-30), 0)
    # V_xc via FFT grid
    gr = to_grid(rho_g, hkl_w)
    rho_r = np.fft.ifftn(gr)*nfft**3   # sum rho_G e^{iGr}
    exc_r, vxc_r = pw92_exc_vxc(rho_r.real)
    vxc_g_full = np.fft.fftn(vxc_r)/nfft**3
    idx = tuple((hkl_w % nfft).T)
    vxc = vxc_g_full[idx]
    E_xc = float(np.sum(rho_r.real*exc_r))*om/nfft**3
    # H and solve
    VHX = np.zeros((ng, ng), complex)
    # potential matrix: V(G-G') lookup needs full conv; reuse FFT: diag in r
    vloc_diag_r = None
    # build local+H+xc in r-space via FFT for matrix-free? -> assemble dense:
    vtot_g = vh + vxc            # on sphere; plus VL matrix for vloc
    gr2 = to_grid(vtot_g, hkl_w)
    vr = np.fft.ifftn(gr2)*nfft**3
    # dense conv matrix from vr: V_{ij} = (1/N^3) sum_r v(r) e^{-i(Gi-Gj)r} =
    # FFT of v at (Gi-Gj)
    vg_full = np.fft.fftn(vr)/nfft**3
    for i in range(ng):
        d = (hkl_w - hkl_w[i]) % nfft
        VHX[:, i] = vg_full[tuple(d.T)]
    Hm = T + VL + VHX + VNL
    w, v = np.linalg.eigh(Hm)
    occ = v[:, :NB]
    # new density on grid
    rho_new_r = np.zeros((nfft,nfft,nfft))
    for bnd in range(NB):
        psi = np.fft.ifftn(to_grid(occ[:,bnd], hkl_w))*nfft**3
        rho_new_r += 2*np.abs(psi)**2/om
    rg_full = np.fft.fftn(rho_new_r)/nfft**3
    rho_new = rg_full[idx]
    drho = np.linalg.norm(rho_new-rho_g)
    rho_g = rho_g + 0.5*(rho_new-rho_g)
    e_band = 2*np.sum(w[:NB])
    e_h = float(0.5*np.real(np.vdot(rho_g, vh)))*om
    e_xci = float(np.real(np.vdot(rho_g, vxc)))*om
    e_loc_v = float(np.real(np.sum(np.conj(rho_g)*np.array([vmap[int(x)] for x in (hkl_w**2).sum(1)]))))*om
    e_one = e_band - 2*e_h - e_xci
    etot = e_one + e_h + E_xc + E_ewald
    print(f"it{it}: eig={w[:4]} e_band={e_band:.5f} E_H={e_h:.5f} E_xc={E_xc:.5f} e_xci={e_xci:.5f} E_loc={e_loc_v:.5f} TOT={etot:.6f} Ha ({etot*2:.5f} Ry) drho={drho:.2e}")
    if drho < 1e-6: break
