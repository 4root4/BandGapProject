#import
import itertools
import numpy as np
import matplotlib.pyplot as plt

# Konstanten
#4H-SiC Gitterkonstanten
A = 3.073
C = 10.053

# Ideale Position von C über Si
U = 3/16

#ħ² durch mₑ
HBAR2_M = 7.619964

#Orbitalenergien in eV von SiC 
ONSITE = {"Si": (-15.050, -5.020, 24.950),
          "C":  (-16.020, -8.970, 23.980)}

# dimensionslose Kopplungen η (V = η·ħ²/(m d²)); ca = s(Si)–p(C), ac = s(C)–p(Si)
ETA = {"ss": -1.57980, "sp_ca": 1.04569, "sp_ac": 1.64281, "pps": 2.70061,
       "ppp": -0.48784, "s*p_ca": 1.53756, "s*p_ac": 4.18428}
NORB = 5                    # s, px, py, pz, s*
VALENCE = {"Si": 4, "C": 4}
CATION = {"Si"}


# Definition der Kristallstruktur und Stapelfolge (4H, abcb) für andere Kristallstrukturen ändern!!!

#abc im Array speichern

SupCell = np.array([[A,0,0], [-A/2, A*np.sqrt(3)/2,0], [0,0,C]])
SI_FRAC = np.array([[0,0,0],[1/3,2/3,1/4], [0,0,1/2],[2/3,1/3,3/4]])

FRAC0 = np.vstack([SI_FRAC,SI_FRAC + [0,0,U]])
SPEC0 = ["Si"] * 4 + ["C"] * 4

#Definiern der Kovalenten Bindungen

def bonds(frac, cell, present, cutoff=2.3):
    out = []
    shifts = np.array(list(itertools.product((-1,0,1), repeat=3)))
    idx = np.where(present)[0]

    for i in idx:
        d = (frac[idx][None] + shifts[:, None] - frac[i]) @ cell
        dist = np.linalg.norm(d, axis=2)
        
        for s, m in zip(*np.where((dist < cutoff) & (dist > 1e-6))):
            out.append((i, idx[m], d[s, m]))
    return out
    

#Definieren des Hamiltonians

def hamilton(si, sj, d):
    r = np.linalg.norm(d)
    l = d / r
    f = HBAR2_M / r**2
    ci = si in CATION   # Atom i ist Si (Kation) oder C (Anion)
    a, b = (ETA["sp_ca"], ETA["sp_ac"]) if ci else (ETA["sp_ac"], ETA["sp_ca"])
    a2, b2 = (ETA["s*p_ca"], ETA["s*p_ac"]) if ci else (ETA["s*p_ac"], ETA["s*p_ca"])
    t = np.zeros((5, 5))
    t[0, 0] = ETA["ss"] * f
    t[0, 1:4], t[1:4, 0] = l * a * f, -l * b * f
    t[1:4, 1:4] = (np.outer(l, l) * (ETA["pps"] - ETA["ppp"]) + np.eye(3) * ETA["ppp"]) * f
    t[4, 1:4], t[1:4, 4] = l * a2 * f, -l * b2 * f
    return t
   

class Hamiltonian:
    def __init__(self, frac, spec, cell, present=None):
        self.present = np.ones(len(spec), bool) if present is None else present
        self.cell = cell
        self.idx = -np.ones(len(spec), int)
        self.idx[self.present] = np.arange(self.present.sum())
        self.n = NORB * self.present.sum()
        self.onsite = np.zeros(self.n)
        for i in np.where(self.present)[0]:
            es, ep, ess = ONSITE[spec[i]]
            self.onsite[self.idx[i] * 5:self.idx[i] * 5 + 5] = (es, ep, ep, ep, ess)
        self.B = bonds(frac, cell, self.present)
        self.T = [hamilton(spec[i], spec[j], d) for i, j, d in self.B]
        self.n_el = sum(VALENCE[s] for s, p in zip(spec, self.present) if p)

    def H(self, k=np.zeros(3)):
        """H(k) mit Bloch-Phasen e^{ik·d}; bei k = 0 reell."""
        H = np.zeros((self.n, self.n), complex if np.any(k) else float)
        for (i, j, d), t in zip(self.B, self.T):
            a, b = self.idx[i] * 5, self.idx[j] * 5
            H[a:a + 5, b:b + 5] += (np.exp(1j * k @ d) if np.any(k) else 1.0) * t
        H[np.diag_indices(self.n)] += self.onsite
        return H

    def recip(self):
        return 2 * np.pi * np.linalg.inv(self.cell).T



# Bänderstruktur eines idealen Kristalles

HSP = {"Γ": [0, 0, 0], "M": [.5, 0, 0], "K": [1/3, 1/3, 0],
       "A": [0, 0, .5], "L": [.5, 0, .5], "H": [1/3, 1/3, .5]}
PATH = ["Γ", "M", "K", "Γ", "A", "L", "H", "A"]


def band_structure(n_seg=60):
    h = Hamiltonian(FRAC0, SPEC0, SupCell)
    B = h.recip()
    xs, Es, ticks, x0 = [], [], [0.0], 0.0
    nseg_tot = len(PATH) - 1
    for iseg, (a, b) in enumerate(zip(PATH[:-1], PATH[1:])):
        ka, kb = np.array(HSP[a]), np.array(HSP[b])
        L = np.linalg.norm((kb - ka) @ B)
        for t in np.linspace(0, 1, n_seg, endpoint=(iseg == nseg_tot - 1)):
            Es.append(np.linalg.eigvalsh(h.H((ka + t * (kb - ka)) @ B)))
            xs.append(x0 + t * L)
        x0 += L
        ticks.append(x0)
    E = np.array(Es)
    nocc = h.n_el // 2
    # Bandkanten zusätzlich auf einem 3D-k-Netz suchen
    vbm, cbm = -np.inf, np.inf
    for k in itertools.product(*(np.arange(m) / m for m in (12, 12, 6))):
        e = np.linalg.eigvalsh(h.H(np.array(k) @ B))
        vbm, cbm = max(vbm, e[nocc - 1]), min(cbm, e[nocc])
    return np.array(xs), E, ticks, nocc, vbm, cbm



# Plotten der Bänder


def plot_bands(x, E, ticks, nocc, vbm, cbm, fn="bandstruktur.png"):
    E = E - vbm
    gap = cbm - vbm
    fig, ax = plt.subplots(figsize=(8, 6.2))
    ax.axhspan(0, gap, color="#f3f2ee")
    for n in range(E.shape[1]):
        ax.plot(x, E[:, n], color="#0072B2" if n < nocc else "#ED3D11", lw=1.3)
    for t in ticks[1:-1]:
        ax.axvline(t, color="#d9d8d3", lw=0.8)
    for y, lab in ((0, ""), (gap, "")):
        ax.axhline(y, color="k", ls="--", lw=0.8)
        ax.text(x[-1], y + 0.15, lab + " ", ha="center", fontsize=8)
    ic = np.argmin(E[:, nocc])
    ax.annotate("", (x[ic], gap), (x[ic], 0), arrowprops=dict(arrowstyle="<->"))
    ax.text(x[ic] + 0.1, gap / 2, f"E$_g$ = {gap:.2f} eV",
            va="center", bbox=dict(fc="white", ec="#d9d8d3"))
    ax.set_xticks(ticks)
    ax.set_xticklabels(PATH)
    ax.set_xlabel("Wellenvektor Γ")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-18, 9)
    ax.set_ylabel("E (eV)")
    ax.set_title("4H-SiC Bänderdiagramm", loc="center")
    fig.savefig(fn, dpi=200, bbox_inches="tight")
    plt.close(fig)




# Ausführen

if __name__ == "__main__":
    x, E, ticks, nocc, vbm, cbm = band_structure()
    print(f"Bandlücke E_g = {cbm - vbm:.3f} eV")
    plot_bands(x, E, ticks, nocc, vbm, cbm)
