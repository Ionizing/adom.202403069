#!/usr/bin/env python3

import numpy as np
# from numpy.typing import NDArray
import matplotlib.pyplot as plt
# from matplotlib.ticker import FuncFormatter
import h5py
# from glob import glob


def lower_triangle_matrix_index(i: int, j: int) -> int:
    if i < j:
        i, j = j, i
    return i * (i + 1) // 2 + j


def get_rtime(isw: int, nsw: int, namdinit: int) -> int:
    return (isw + namdinit) % (nsw - 2)


class AveragedResult:
    """
    DATASET "delta_et",                         [npairs, nsw]
    DATASET "eigs_t",                           [namdtime, nbasis]
    DATASET "ndigit",
    DATASET "phonon_spectra_frequencies",       [nfreq]
    DATASET "phonons_absp_t",                   [namdtime, nfreq]
    DATASET "phonons_emit_t",                   [namdtime, nfreq]
    DATASET "phonons_spectra",                  [npairs, nfreq]
    DATASET "photon_spectra_xvals",             [npoints]
    DATASET "photons_absp_t",                   [namdtime, npoints]
    DATASET "photons_emit_t",                   [namdtime, npoints]
    DATASET "potim",
    DATASET "proj_t",                           [namdtime, nbasis, nions, nproj]
    DATASET "prop_energy",                      [namdtime]
    DATASET "psi_t",                            [namdtime, nbasis]
    DATASET "sh_phonons_t",                     [namdtime, nbasis, nbasis]
    DATASET "sh_photons_t",                     [namdtime, nbasis, nbasis]
    DATASET "sh_energy",                        [namdtime]
    DATASET "sh_pops",                          [namdtime, nbasis]
    DATASET "time",                             [namdtime]

    where
        npairs  = nbasis*(nbasis+1) / 2
        npoints = points sampled in the photon frequency domain
        nproj = 9 (by default)
    """
    def __init__(self, fname:str="averaged_results.h5"):
        with  h5py.File("averaged_results.h5") as f:
            for k in f:
                setattr(self, k, f[k][()])
        self.basis_labels = self.basis_labels.tobytes().decode("utf-8").split()
        np.savez("shdata.npz",
                 eigs_t=self.eigs_t,
                 proj_t=np.sum(self.proj_t, axis=(2,3)),
                 psi_t=self.psi_t,
                 prop_energy=self.prop_energy,
                 sh_energy=self.sh_energy,
                 sh_pops=self.sh_pops,
                 time=self.time)
        pass


    def plot_wfn_propagation(self, _ax=None):
        """
        TODO
        """
        if _ax is None:
            fig, ax = plt.subplots(figsize=(6,4))

        namdtime, nbasis = self.psi_t.shape
        prop_energy = self.prop_energy
        eigs_t = self.eigs_t
        time = self.time
        T = np.array([time for _ in range(nbasis)])
        c = self.psi_t * np.sum(self.proj_t, axis=(2, 3))
        kmap = ax.scatter(T.T, eigs_t, c=c,
                          cmap="Reds", s=15, lw=0.0,
                          rasterized=True, vmin=0, vmax=1)
        # ax.plot(T[0], prop_energy)
        ax.set_xlim(time.min(), time.max())
        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("|C(t)|²")
        ax.set_title("wfn propagation")

        if _ax is None:
            cb = fig.colorbar(kmap, fraction=0.046, pad=0.01)
            fig.tight_layout(pad=0.5)
            print("Writing wfn_propagation.png ...")
            fig.savefig("wfn_propagation.png", dpi=400)
        pass

    
    def plot_surfhop(self, _ax=None):
        """
        TODO
        """

        if _ax is None:
            fig, ax = plt.subplots(figsize=(6,4))
        
        namdtime, nbasis = self.psi_t.shape
        sh_energy = self.sh_energy
        eigs_t = self.eigs_t
        time = self.time
        T = np.array([time for _ in range(nbasis)])
        c = self.sh_pops * np.sum(self.proj_t, axis=(2, 3))
        kmap = ax.scatter(T.T, eigs_t, c=c,
                          cmap="Reds", s=15, lw=0.0,
                          rasterized=True, vmin=0, vmax=1)
        # ax.plot(T[0], sh_energy)
        ax.set_xlim(time.min(), time.max())
        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("E - Ef (eV)")
        ax.set_title("surface hopping population")

        if _ax is None:
            cb = fig.colorbar(kmap, fraction=0.046, pad=0.01)
            fig.tight_layout(pad=0.5)
            print("Writing surfhop.png ...")
            fig.savefig("surfhop.png", dpi=400)
        pass


    def plot_phonon_spectra(self, pair, _ax=None):
        """
        Plot the phonon spectra for given traisition denoted by the band pair i -> j.
        The diagonal part i -> i is set to be the band eigenvalue itself.

        The band pair should count from 0 and match the indices of current basis space.
        """
        freq = self.phonon_spectra_frequencies
        spectra = self.phonons_spectra

        if _ax is None:
            fig, ax = plt.subplots(figsize=(6,4))
        else:
            ax = _ax

        idx = lower_triangle_matrix_index(*pair)
        ax.plot(freq, spectra[idx])
        ax.set_xlim(0, 3000)
        ax.set_xlabel("Wavenumber (cm-1)")
        ax.set_ylabel("Intensity (arb. unit.)")
        ax.set_title("Phonon spectra")

        if _ax is None:
            fig.tight_layout(pad=0.5)
            print("Writing phonon_spectra.png ...")
            fig.savefig("phonon_spectra.png", dpi=400)
        pass


    def plot_phonon_waterfall(self, times, _ax=None):
        """
        Plot waterfall diagram of time-resolved phonon spectra for all transitions.
        """
        potim = self.potim
        freq = self.phonon_spectra_frequencies
        emitted = self.phonons_emit_t
        absorped = self.phonons_absp_t

        if _ax is None:
            fig, ax = plt.subplots(figsize=(6,4))
        else:
            ax = _ax

        for (i, t) in enumerate(times):
            ax.plot(freq, emitted[t,:] + i/5, label="phon. emit. t={}ps".format(t*potim/1000))
            ax.plot(freq, absorped[t,:] - i/5, label="phon. absp. t={}ps".format(t*potim/1000))

        ax.set_yticks([])
        ax.set_xlim(0, 3000)
        # ax.legend()
        ax.set_xlabel("Wavenumber (cm-1)")
        ax.set_ylabel("Intensity (arb. unit.)")
        ax.set_title("Time-resolved phonon spectra")
        if _ax is None:
            fig.tight_layout(pad=0.5)
            print("Writing phonon_waterfall.png ...")
            fig.savefig("phonon_waterfall.png", dpi=400)
        pass


    def plot_phonon_photon_contribution(self,
                                        *,
                                        iframes=None,
                                        pngfname="contributions.png",
                                        _ax=None):
        sh_phonons_t = self.sh_phonons_t * 1E5
        sh_photons_t = self.sh_photons_t * 1E5
        assert sh_phonons_t.shape == sh_photons_t.shape

        namdtime, nbasis = self.psi_t.shape
        contributions = np.zeros_like(sh_phonons_t)

        for i in range(nbasis):
            for j in range(i+1, nbasis):
                # upper triangle matrix is phonon contribution
                contributions[:, i, j] = sh_phonons_t[:,i,j] + sh_phonons_t[:,j,i]

                # lower triangle matrix is photon contribution
                contributions[:, j, i] = sh_photons_t[:,i,j] + sh_photons_t[:,j,i]

        if iframes is None:
            contributions = np.mean(contributions, axis=0)
        else:
            contributions = np.mean(contributions[iframes, :, :], axis=0)

        if _ax is None:
            fig, ax = plt.subplots(figsize=(5, 4))
        else:
            ax = _ax
        
        contributions_abs = np.abs(contributions)
        valrange = contributions_abs.max() - contributions_abs.min()
        # vmin = contributions.min() - valrange * 0.15
        vmax = contributions_abs.max() + valrange * 0.2
        img = ax.pcolormesh(contributions_abs, cmap="Reds",
                            linewidth=0,
                            aa=True,
                            # vmin=vmin,
                            vmax=vmax,
                            edgecolor="none")

        if nbasis <= 15:
            for (i, j), z in np.ndenumerate(contributions):
                ax.text(j+0.5, i+0.5, "{:0.2f}".format(z),
                        ha="center", va="center",
                        color="k")

        if _ax is None:
            cb = fig.colorbar(img, fraction=0.046, pad=0.01)
            cb.ax.set_title("x10$^{-5}$")

            ticks = [i+0.5 for i in range(nbasis)]
            # ticklabels = [
                    # "vK-d", "vK+u", "vK+d", "vK-u",
                    # "cK+d", "cK-u", "cK-d", "cK+u",
                    # ]
            ticklabels = self.basis_labels
            ax.set_xticks(ticks)
            ax.set_xticklabels(ticklabels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(ticklabels)

            fig.tight_layout(pad=0.5)
            fig.suptitle("Photon(↖) and phonon(↘) contributions.")
            print(f"Writing {pngfname} ...")
            fig.savefig(pngfname, dpi=400)
        
        return contributions


    def plot_photon_waterfall(self, times, _ax=None):
        """
        Plot waterfall diagram of time-resolved photon spectra for all transitions.
        """
        potim = self.potim
        freq = self.photon_spectra_xvals
        emitted = self.photons_emit_t
        absorped = self.photons_absp_t

        if _ax is None:
            fig, ax = plt.subplots(figsize=(6,4))
        else:
            ax = _ax

        for (i, t) in enumerate(times):
            ax.plot(freq, emitted[t,:], label="phot. emit. t={}ps".format(t*potim/1000))
            ax.plot(freq, absorped[t,:], label="phot. absp. t={}ps".format(t*potim/1000))

        ax.set_xlim(0, freq.max())
        ax.set_yticks([])
        ax.set_xlabel("Photon energy (eV)")
        ax.set_ylabel("Intensity (arb. unit.)")
        ax.set_title("Time-resolved photon spectra")
        if _ax is None:
            fig.tight_layout(pad=0.5)
            print("Writing photon_waterfall.png ...")
            fig.savefig("photon_waterfall.png", dpi=400)
        pass


    def plot_valley_population(self, pngfname="valley_populations.png"):
        # valleys_with_spin = ["vK-d", "vK+u", "vK+d", "vK-u",
                             # "cK+d", "cK-u", "cK-d", "cK+u",]
        valleys_with_spin = self.basis_labels
        namdtime, nbasis = self.psi_t.shape
        time = self.time
        pops = self.sh_pops


        fig, axs = plt.subplots(nrows=2, figsize=(5, 4), height_ratios=[1, 0.6], sharex=True)
        axs[0].plot(time, pops[:,4:], ls="-",  label=valleys_with_spin[4:])
        axs[0].plot(time, pops[:,:4], ls="--", label=valleys_with_spin[:4])
        axs[0].set_xlabel("Time (fs)")
        axs[0].set_ylabel("Population")
        axs[0].set_title("Population for each band (with spin)")
        axs[0].legend()

        pops_valley = np.sum(pops[:,[[0,0],[1,4],[2,3]]], axis=-1)
        pops_valley[:,0] /= 2
        # valleys_cbvb_only = [valleys_with_spin[0], "CB@K+", "CB@K-"]

        polarization = (pops[:,5] - pops[:,4]) / (pops[:,4] + pops[:,5] + 2.5E-5)
        axs[1].plot(time, polarization, label="electron", color="red")
        axs[1].axhline(y=0, ls='--', lw=0.5)

        hpops = 1 - pops
        hpolarization = (hpops[:,3] - hpops[:,2]) / (hpops[:,3] + hpops[:,2] + 2.5E-5)
        axs[1].plot(time, hpolarization, label="hole", color="blue")

        # axs[1].plot(time, pops_valley[:,1:], label=valleys_cbvb_only[1:])
        # axs[1].plot(time, pops_valley[:,0 ], ls="--", label=valleys_cbvb_only[0])
        axs[1].set_xlabel("Time (fs)")
        axs[1].set_ylabel("Polarization")
        axs[1].set_ylim(-1.05, 1.05)
        axs[1].legend()

        fig.tight_layout(pad=0.5)
        print(f"Writing {pngfname} ...")
        fig.savefig(pngfname, dpi=400)
        pass


if '__main__' == __name__:
    ar = AveragedResult()
    ar.plot_wfn_propagation()
    ar.plot_surfhop()
    # ar.plot_phonon_spectra(pair=(5,6))
    # ar.plot_phonon_waterfall(list(range(1000, 5000, 1000)))
    # ar.plot_photon_waterfall(list(range(1000, 5000, 1000)))
    ar.plot_valley_population()
    # ar.plot_phonon_photon_contribution(
        # pngfname="all_contributions.png"
    # )
    # ar.plot_phonon_photon_contribution(
        # iframes=slice(0, 200),
        # pngfname="first_200fs_contributions.png"
    # )
    # ar.plot_phonon_photon_contribution(
        # iframes=slice(0, 400),
        # pngfname="first_400fs_contributions.png"
    # )
    # ar.plot_phonon_photon_contribution(
        # iframes=slice(1000, 5000),
        # pngfname="first_4000fs_contributions.png"
    # )
    pass
