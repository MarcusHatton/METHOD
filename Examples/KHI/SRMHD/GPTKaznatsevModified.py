import numpy as np
import matplotlib.pyplot as plt
import h5py
import pickle

def mhd_spherical_energy_spectra(frame, bins=50):
    vx = frame['Primitive/vx'][:]
    vy = frame['Primitive/vy'][:]
    vz = frame['Primitive/vz'][:]
    Bx = frame['Primitive/Bx'][:]
    By = frame['Primitive/By'][:]
    Bz = frame['Primitive/Bz'][:]
    rho = frame['Primitive/rho'][:]
    W = frame['Auxiliary/W'][:]


    N = vx.shape[0]
    assert all(f.shape == (N, N, N) for f in [vy, vz, Bx, By, Bz]), "Fields must be same 3D shape"

    # FFT of each component
    def fft3(f):
        return np.fft.fftshift(np.fft.fftn(f))

    rhof, Wf = fft3(rho), fft3(W)
    vxf, vyf, vzf = fft3(vx), fft3(vy), fft3(vz)
    Bxf, Byf, Bzf = fft3(Bx), fft3(By), fft3(Bz)

    # Compute energy density per mode
    # E_kin_density = 0.5 * (np.abs(vxf)**2 + np.abs(vyf)**2 + np.abs(vzf)**2)
    E_kin_density = rhof * Wf
    E_mag_density = 0.5 * (np.abs(Bxf)**2 + np.abs(Byf)**2 + np.abs(Bzf)**2)

    # Create k-space grid
    k = np.fft.fftshift(np.fft.fftfreq(N)) * N
    KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij')
    k_mag = np.sqrt(KX**2 + KY**2 + KZ**2)

    # Flatten
    k_flat = k_mag.flatten()
    E_kin_flat = E_kin_density.flatten()
    E_mag_flat = E_mag_density.flatten()

    # Bin by k magnitude
    k_bins = np.linspace(0, k_mag.max(), bins + 1)
    k_vals = 0.5 * (k_bins[:-1] + k_bins[1:])
    E_kin_1D = np.zeros(bins)
    E_mag_1D = np.zeros(bins)

    bin_indices = np.digitize(k_flat, k_bins) - 1
    for i in range(bins):
        mask = bin_indices == i
        if np.any(mask):
            E_kin_1D[i] = E_kin_flat[mask].mean()
            E_mag_1D[i] = E_mag_flat[mask].mean()

    return k_vals, E_kin_1D, E_mag_1D

if __name__ == "__main__":
    fs = []
    n_files = 8
    for n in range(n_files):
        fs.append(h5py.File(f'./3d/KHI/dp_200x200x200_{n}.hdf5', 'r'))
    # print(fs)
    nx = ny = nz = 200

    KE_Spec = 0
    n_slice = 7
    frame = fs[n_slice]
    # Compute spectra
    k_vals, E_kin, E_mag = mhd_spherical_energy_spectra(frame)

    with open('./pickles/KEMESpec3D_200x200x200_t=20.pickle', 'wb') as filehandle:
        pickle.dump([k_vals, E_kin, E_mag], filehandle)

    # Plot
    plt.figure(figsize=(8,6))
    plt.loglog(k_vals, E_kin, label="Kinetic Energy", lw=2)
    plt.loglog(k_vals, E_mag, label="Magnetic Energy", lw=2)

    # Overlay power-law lines (Kolmogorov, Iroshnikov-Kraichnan)
    valid = k_vals > 0
    k_sample = k_vals[valid]
    k0 = k_sample[len(k_sample)//2]
    E0 = E_kin[valid][len(k_sample)//2]
    plt.loglog(k_sample, E0*(k_sample/k0)**(-5/3), 'k--', label=r"$k^{-5/3}$")
    plt.loglog(k_sample, E0*(k_sample/k0)**(-3/2), 'k-.', label=r"$k^{-3/2}$ (IK)")

    # Final touches
    plt.xlabel("Wavenumber $k$")
    plt.ylabel("Energy Spectrum")
    plt.title("MHD Energy Spectra")
    plt.legend()
    plt.grid(True, which="both", ls=":")
    plt.tight_layout()
    plt.savefig(f'./plots/KEMESpec3D_200x200x200_t=20.pdf', bbox_inches='tight')
    plt.show()
