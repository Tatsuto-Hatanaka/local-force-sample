import os
from time import time
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# print(np.__version__) 1.26.1
# print(sp.__version__) 1.11.3

class local_force_simple_model:
    def __init__(self, beta=100, D=1, U=2, e_num=501, m_init=0.1, m_min=1e-8, m_max=1.5):
        if not os.path.isdir("figures"): os.mkdir("figures")
        self.D  = D
        self.U  = U
        e_max = D + U*m_max/2 + 0.5
        self.e_range = np.linspace(-e_max, e_max, e_num) + 1e-10 + 0j
        # for the self-consistent loop
        self.mag_init = 1
        self.mag_min = m_min
        self.mag_max = m_max
        mag = 1.0
        gf0   = self.calc_green_function(self.e_range)
        gf_up = self.calc_green_function(self.e_range + self.U * mag/2)
        gf_dn = self.calc_green_function(self.e_range - self.U * mag/2)
        dos0   = - gf0.imag/np.pi
        dos_up = - gf_up.imag/np.pi
        dos_dn = - gf_dn.imag/np.pi
        integrated_dos0    = self.cumulative_simpson(dos0, self.e_range.real)
        integrated_dos_up  = self.cumulative_simpson(dos_up, self.e_range.real)
        integrated_dos_dn  = self.cumulative_simpson(dos_dn, self.e_range.real)
        integrated_dos_all = self.cumulative_simpson(dos_up + dos_dn, self.e_range.real)
        fig, ax = plt.subplots(2,2, figsize=(12,8), tight_layout=True)
        ax[0,0].set_title("DOS"+r"($\epsilon$)")
        ax[0,0].set_xlabel(r"$\epsilon$")
        ax[0,0].plot(self.e_range.real, dos0, c="black", label=r"$D^{0}$")
        ax[0,0].plot(self.e_range.real, dos_up, c="blue", label=r"$D^{\uparrow}$")
        ax[0,0].plot(self.e_range.real, dos_dn, c="red", label=r"$D^{\downarrow}$")
        ax[0,1].set_title("N"+r"($\epsilon$)")
        ax[0,1].set_xlabel(r"$\epsilon$")
        ax[0,1].plot(self.e_range.real, integrated_dos0, c="black", label=r"$N^{0}$")
        ax[0,1].plot(self.e_range.real, integrated_dos_up, c="blue", label=r"$N^{\uparrow}$")
        ax[0,1].plot(self.e_range.real, integrated_dos_dn, c="red", label=r"$N^{\downarrow}$")
        ax[1,0].set_title("Re[G"+r"($\epsilon$)"+"]")
        ax[1,0].set_xlabel(r"$\epsilon$")
        ax[1,0].plot(self.e_range.real, gf0.real, c="black", label=r"Re[$G^{0}$]")
        ax[1,0].plot(self.e_range.real, gf_up.real, c="blue", label=r"Re[$G^{\uparrow}$]")
        ax[1,0].plot(self.e_range.real, gf_dn.real, c="red", label=r"Re[$G^{\downarrow}$]")
        ax[1,1].set_title(r"Im[$G(\epsilon)$]")
        ax[1,1].set_xlabel(r"$\epsilon$")
        ax[1,1].plot(self.e_range.real, gf0.imag, c="black", label=r"Im[$G^{0}$]")
        ax[1,1].plot(self.e_range.real, gf_up.imag, c="blue", label=r"Im[$G^{\uparrow}$]")
        ax[1,1].plot(self.e_range.real, gf_dn.imag, c="red", label=r"Im[$G^{\downarrow}$]")
        text = r"$U=$"+"{}\n".format(self.U) + r"$m_{0}=$"+"{}".format(mag)
        for ia1,ia2 in [[ia1, ia2] for ia1 in range(2) for ia2 in range(2)]:
            ax[ia1,ia2].axvline(0, color="black", linestyle='dashed', linewidth=0.5)
            ax[ia1,ia2].axhline(0, color="black", linestyle='dashed', linewidth=0.5)
            ax[ia1,ia2].set_xlim(-e_max, e_max)
            if ia1==0: ax[ia1,ia2].set_ylim(-1.0, 1.0)
            elif ia1==1: ax[ia1,ia2].set_ylim(-3.0, 3.0)
            ax[ia1,ia2].text(x=0.85, y=0.05, s=text, c="black", fontsize=13, ha="center"\
                        , transform=ax[ia1,ia2].transAxes, in_layout=True)
            ax[ia1,ia2].legend()
        fig.savefig("figures/dos_n_gfs.pdf")

        self.mags = np.zeros(e_num)
        self.j0s  = np.zeros(e_num)
        for ie, e in enumerate(self.e_range):
            fermi_dist = 0.5 * (1 - np.tanh(0.5 * beta * (self.e_range - e)))
            self.mags[ie] = self.magnetization_self_consistent(ie, fermi_dist)
            self.gf_up = self.calc_green_function(self.e_range + self.U * self.mags[ie]/2)
            self.gf_dn = self.calc_green_function(self.e_range - self.U * self.mags[ie]/2)
            self.j0s[ie] = self.calc_j0(ie, fermi_dist)

    def calc_green_function(self, e_range):
        gfs = (2/self.D) * (e_range/self.D - np.sign(e_range) * np.sqrt((e_range/self.D)**2-1))
        return gfs

    def magnetization_self_consistent(self, ie, fermi_dist):
        try:
            solver = sp.optimize.root_scalar(f=self.diff_magnetization, args=(ie, fermi_dist), method="brentq"\
                                , bracket=(self.mag_min,self.mag_max), x0=self.mag_init, xtol=1e-4, rtol=1e-4, maxiter=1000)
            mag = solver.root
        except:
            mag = 0.0
        return mag

    def diff_magnetization(self, mag, ie, fermi_dist):
        gf_up = self.calc_green_function(self.e_range + self.U * mag/2)
        gf_dn = self.calc_green_function(self.e_range - self.U * mag/2)
        mag_from_gf = - sp.integrate.simpson(fermi_dist.real*(gf_up-gf_dn).imag, self.e_range.real)/np.pi
        return mag - mag_from_gf

    def calc_j0(self, ie, fermi_dist):
        integrand = fermi_dist.real * ((self.U*self.mags[ie])*(self.gf_up-self.gf_dn)\
                    + (self.U*self.mags[ie])**2 * self.gf_up * self.gf_dn)
        j0 = - sp.integrate.simpson(integrand.imag, self.e_range.real)/(4*np.pi)
        return j0

    def cumulative_simpson(self, vals, e_range):
        assert len(vals)==len(e_range)
        n     = len(e_range)
        start = 0
        stop  = n-2 if n%2==1 else n-3
        de    = (e_range[-1]-e_range[0])/(n-1)
        outs  = np.zeros(n)
        outs[start:stop:2]     += vals[start:stop:2]
        outs[start+1:stop+1:2] += vals[start+1:stop+1:2]
        outs[start+2:stop+2:2] += vals[start+2:stop+2:2]
        outs *= de/3.
        # https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_rule_for_irregularly_spaced_data
        if n%2==0:
            alpha = 5/12
            beta  = 2/3
            eta   = 1/12
            outs[-1] += alpha * vals[-1]*de
            outs[-2] += beta * vals[-2]*de
            outs[-3] -= eta * vals[-3]*de
        return np.cumsum(outs)


if __name__=="__main__":
    l = local_force_simple_model(U=2.0)
    fig, ax = plt.subplots(1,2, figsize=(12,5), tight_layout=True)
    ax[0].set_title(r"$m_0(\mu)$")
    ax[0].set_xlabel(r"$\mu$")
    ax[0].plot(l.e_range.real, l.mags, label="magnetization")
    ax[1].set_title(r"$J_0(\mu)$")
    ax[1].set_xlabel(r"$\mu$")
    ax[1].plot(l.e_range.real, l.j0s, label=r"$J_0$")
    for ia in range(2):
        ax[ia].set_xlim(-l.U, l.U)
        ax[ia].axvline(0, color="black", linestyle='dashed', linewidth=0.5)
        ax[ia].axhline(0, color="black", linestyle='dashed', linewidth=0.5)
        ax[ia].legend()
    fig.savefig("figures/mag_j0.pdf")