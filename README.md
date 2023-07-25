# BSM soil model python package

This package contains the python implementation of the BSM soil model. The model is described in the following paper:

[Hyperspectral radiative transfer modeling to explore the combined retrieval of biophysical parameters and canopy fluorescence from FLEX – Sentinel-3 tandem mission multi-sensor data](https://www.sciencedirect.com/science/article/pii/S0034425717303607)

## Installation
```bash
pip install https://github.com/MarcYin/BSM_soil/archive/refs/heads/main.zip
```

The model is implemented in the `bsm` module. The `bsm` module contains the following parameters and their ranges:

| Parameter | Description | Range |
| --- | --- | --- |
| B | Soil brightness | 0.0-1.0 |
| lat | Latitude | 20°–40° |
| lon | Longitude | 45°–65° |
| SMp | Soil moisture | 0-55 |

Example usage:
```python
from BSM import BSM

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
```
