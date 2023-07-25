import os
import numpy as np
from typing import Tuple

# Defining global constants and lookup tables.
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

LOOKUP_TABLE = np.array([
    1, 1, 2, 6, 24, 120, 720, 5040, 40320,
    362880, 3628800, 39916800, 479001600,
    6227020800, 87178291200, 1307674368000,
    20922789888000, 355687428096000, 6402373705728000,
    121645100408832000, 2432902008176640000], dtype='int64')


def load_BSM_paras() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load BSM parameters from a data file.

    :return: A tuple containing GSV, nw, kw numpy arrays.
    """
    data = np.loadtxt(os.path.join(data_dir, 'BSM_soil_paras.txt'), skiprows=1, dtype=float, delimiter=',')
    return data[:, 0:3], data[:, [3]], data[:, [4]]


BSM_paras = load_BSM_paras()


def fast_factorial(n: int) -> int:
    """
    Calculates the factorial of a number using lookup table.

    :param n: The number to find the factorial of.
    :return: The factorial of n.
    """
    if n > 20:
        raise ValueError
    return LOOKUP_TABLE[n]

def tav(alpha: float, nr: float) -> float:
    """
    Computes the reflectance for alpha degrees using Stern's formula (Lekner & Dorf, 1988).
    
    :param alpha: The alpha degrees.
    :param nr: The refraction index.
    :return: The computed reflectance value.
    """
    n2 = nr**2
    np_ = n2 + 1  
    nm_ = n2 - 1 
    a = (nr + 1)**2 / 2
    k = -(n2 - 1)**2 / 4

    sin_a = np.sin(np.deg2rad(alpha))
    if alpha != 0:
        B2 = sin_a**2 - np_ / 2
        if alpha == 90:
            B1 = B2 * 0
        else:
            B1 = np.sqrt(B2**2 + k)
        b = B1 - B2
        b3 = b**3
        a3 = a**3
        ts = (k**2 / (6 * b3) + k / b - b / 2) - (k**2 / (6 * a3) + k / a - a / 2)
        tp1 = -2 * n2 * (b - a) / np_**2
        tp2 = -2 * n2 * np_ * np.log(b / a) / nm_**2
        tp3 = n2 * (1 / b - 1 / a) / 2
        tp4 = 16 * n2**2 * (n2**2 + 1) * np.log((2 * np_ * b - nm_**2) / (2 * np_*a - nm_**2)) / (np_**3 * nm_**2)
        tp5 = 16 * n2**2 * n2 * (1 / (2 * np_ * b - nm_**2) - 1 / (2 * np_ * a - nm_**2)) / np_**3
        tp = tp1 + tp2 + tp3 + tp4 + tp5
        Tav = (ts + tp) / (2 * sin_a**2)
    else:
        Tav = 4 * nr / (nr + 1)**2
    return Tav

def soilwat(rdry: np.ndarray, nw: np.ndarray, kw: np.ndarray, SMp: float, SMC: float, deleff: float) -> np.ndarray:
    """
    Calculates the wet reflectance based on the dry reflectance, refractive index of water,
    absorption coefficient of water, soil moisture, empirical parameter for BSM, and film thickness.
    from https://github.com/Christiaanvandertol/SCOPE/blob/master/src/RTMs/BSM.m

    :param rdry: Dry reflectance values.
    :param nw: Refractive index of water.
    :param kw: Absorption coefficient of water.
    :param SMp: Soil moisture parameter.
    :param SMC: Soil moisture content, an empirical parameter for BSM.
    :param deleff: Effective film thickness.
    :return: Wet reflectance.
    """
    k = np.arange(7)
    nk = len(k)
    mu = (SMp - 5) / SMC
    if mu <= 0:
        rwet = rdry
    else:
        # Lekner & Dorf (1988) modified soil background reflectance
        rbac = 1 - (1-rdry) * (rdry * tav(90, 2.0 / nw) / tav(90, 2.0) + 1 - rdry) # Rbac
        p = 1 - tav(90, nw) / nw**2   # rho21, water to air, diffuse

        # reflectance of water film top surface, use 40 degrees incidence angle,
        # like in PROSPECT
        Rw = 1 - tav(40, nw)  # rho12, air to water, direct

        fmul = np.zeros((nk, 1))
        for i in range(len(k)):
            ret = np.exp(-mu) * mu**i / fast_factorial(i)
            fmul[i] = ret
        tw = np.exp(-2 * kw * deleff * k)  # two-way transmittance, exp(-2*kw*k Delta)
        Rwet_k = Rw + (1 - Rw) * (1 - p) * tw * rbac / (1 - p * tw * rbac)
        rwet = rdry * fmul[0] + Rwet_k[:, 1:].dot(fmul[1:])
    return rwet


def BSM(B: float, lat: float, lon: float, SMp: float) -> np.ndarray:
    """
    This function computes the Brightness Soil Moisture (BSM) from input parameters.

    :param B: Brightness parameter
    :param lat: Latitude in degrees
    :param lon: Longitude in degrees
    :param SMp: Soil moisture parameter
    :return: Wet reflectance (rwet) as a numpy array
    """
    SMC = 25  # empirical parameter (fixed) for BSM
    film = 0.015  # empirical parameter (fixed) for BMS
    GSV, nw, kw = BSM_paras

    latRad = np.deg2rad(lat)
    lonRad = np.deg2rad(lon)
    f1 = B * np.sin(latRad)
    f2 = B * np.cos(latRad) * np.sin(lonRad)
    f3 = B * np.cos(latRad) * np.cos(lonRad)
    rdry = f1 * GSV[:, 0] + f2 * GSV[:, 1] + f3 * GSV[:, 2]
    rdry = np.expand_dims(rdry, axis=1)

    rwet = soilwat(rdry, nw, kw, SMp, SMC, film)

    return rwet

if __name__ == '__main__':
    import time
    B = 0.5 # 0.-1.0
    lat = 30.0 # 20°–40°
    lon = 80.0 # 45°–65°
    SMp = 5.0 # 0-55
    start = time.time()
    rwet = BSM(B, lat, lon, SMp)
    print(f'Elapsed time: {time.time() - start}')
    import matplotlib.pyplot as plt
    plt.figure()
    plt.title(f'B={B}, lat={lat}, lon={lon}, SMp={SMp}')
    plt.plot(range(400, 2401,), rwet)
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Reflectance')
    plt.show()