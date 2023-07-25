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