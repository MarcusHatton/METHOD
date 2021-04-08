def combo():
    """
    Create sum of sine waves called total
    """
    return totalplot, total

plot, waveform = combo()

solution = np.fft.rfft(waveform)

matplotlib.plot(solution)