import numpy as np
import matplotlib.pyplot as plt
import h5py
import pickle

def spherical_average_1d_fft(data, bins=50):
    N = data.shape[0]
    assert data.shape == (N, N, N), "Data must be cubic 3D array"

    # Compute the 3D Fourier transform
    fft_data = np.fft.fftn(data)
    fft_data = np.fft.fftshift(fft_data)  # Center zero frequency

    # Compute power spectrum (magnitude squared)
    power = np.abs(fft_data)**2

    # Create grid of frequencies
    kx = np.fft.fftshift(np.fft.fftfreq(N)) * N
    ky = np.fft.fftshift(np.fft.fftfreq(N)) * N
    kz = np.fft.fftshift(np.fft.fftfreq(N)) * N
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    k_mag = np.sqrt(KX**2 + KY**2 + KZ**2)

    # Flatten the arrays
    k_flat = k_mag.flatten()
    p_flat = power.flatten()

    # Bin by k magnitude
    k_bins = np.linspace(0, k_mag.max(), bins + 1)
    k_vals = 0.5 * (k_bins[1:] + k_bins[:-1])
    power_spectrum = np.zeros(bins)
    counts = np.zeros(bins)

    # Bin the data
    bin_indices = np.digitize(k_flat, k_bins) - 1
    for i in range(bins):
        mask = bin_indices == i
        if np.any(mask):
            power_spectrum[i] = p_flat[mask].mean()
            counts[i] = mask.sum()

    return k_vals, power_spectrum

def GetKESF(frame):
    """
    Retrieves and computes the kinetic energy density for each frame in a single fluid animation.
    Parameters
    ----------
    anim : object
        animation class containing all user def variables
    frame : Array
        Frame from the animation class containing all user def variables at the time we want
    """
    rho = frame['Primitive/rho'][:]
    W = frame['Auxiliary/W'][:]
    KE = rho * W * (W-1)
    return KE

if __name__ == "__main__":
    fs = []
    n_files = 11
    for n in range(n_files):
        fs.append(h5py.File(f'./3d/KHI/dp_80x80x80_{n}.hdf5', 'r'))
    # print(fs)

    nx = ny = nz = 80
    
    KE_Spec = 0
    n_slice = 2
    frame = fs[n_slice]
    KE_Spec = GetKESF(frame)
    k_vals, Pk = spherical_average_1d_fft(KE_Spec)

    with open(f'./pickles/KESpec_{nx}x{ny}x{nz}_t{n_slice*3}.pickle', 'wb') as filehandle:
        pickle.dump([k_vals, Pk], filehandle)

    # Plot power spectrum
    plt.figure(figsize=(8, 6))
    plt.loglog(k_vals, Pk, label="Measured $P(k)$", lw=2)

    # Overlay a power-law line ∝ k^(-5/3)
    # Normalize the power-law to match the magnitude of Pk in mid-k range
    valid = (k_vals > 0)
    k_sample = k_vals[valid]
    Pk_sample = Pk[valid]

    # Choose a region (e.g., middle 1/3 of nonzero k's) to normalize
    mid_idx = len(k_sample) // 2
    k0 = k_sample[mid_idx]
    Pk0 = Pk_sample[mid_idx]
    power_law = Pk0 * (k_sample / k0)**(-5/3)

    plt.loglog(k_sample, power_law, 'k--', label=r"$k^{-5/3}$", lw=2)

    # Final touches
    plt.xlabel("Wavenumber $k$")
    plt.ylabel("Power Spectrum $P(k)$")
    plt.title("1D Spherical Power Spectrum with $k^{-5/3}$ Scaling")
    plt.legend()
    plt.grid(True, which="both", ls=":")
    plt.tight_layout()
    plt.savefig(f'plots/KESpec3D_{nx}x{ny}x{nz}_t{n_slice*3}.pdf', bbox_inches='tight')
    plt.show()
