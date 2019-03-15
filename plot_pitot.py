#!/usr/bin/env python3

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np

def pitot_velocity(dp, rho):
    """Given pressure difference dp [Pa] and air density
    rho [kg/m^3], returns the wind speed measured by the
    pitot tube."""
    return np.sqrt(2 * dp / rho)

def running_mean(x, n):
    """Running mean with the window n."""
    return np.convolve(x, np.ones((n,)) / n, mode='same')

def read_pitot_from_toa5(filenames):
    if type(filenames) is str:
        print('Reading ', filenames)
        data = [line.rstrip() for line in open(filenames).readlines()[4:]]
    elif type(filenames) is list:
        data = []
        for filename in filenames:
            print('Reading ', os.path.basename(filename))
            data += [line.rstrip() for line in open(filename).readlines()[4:]]
    else:
        raise RuntimeError('filenames must be string or list')

    dp, times = [], []
    for line in data:
        line = line.replace('"', '').split(',')
        t = line[0]
        if len(t) == 19:
            time = datetime.strptime(t, '%Y-%m-%d %H:%M:%S')
        else:
            time = datetime.strptime(t[:19], '%Y-%m-%d %H:%M:%S')
            time += timedelta(seconds=float(t[-2:]))
        times.append(time)
        dp.append(float(line[2]))
    return np.array(times), np.array(dp) * 1e3


time, dp = read_pitot_from_toa5('TOA5_SUSTAINpresMKSX4.pressure_8_2019_03_14_1900.dat')

t0 = datetime(2019, 3, 14, 19, 4)
t1 = datetime(2019, 3, 14, 19, 6)

mask = (time >= t0) & (time <= t1)

time, dp = time[mask], dp[mask]

dp_offset = 0.5 * (np.mean(dp[:10]) + np.mean(dp[-10:]))
dp -= dp_offset

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111, xlim=(t0, t1), ylim=(-60, 80))
plt.plot(time, dp, lw=0.4)
plt.plot(time, running_mean(dp, 50), 'r-', lw=1)
plt.plot([t0, t1], [0, 0], 'k--')
plt.grid()
plt.xlabel('Time')
plt.ylabel('Pressure difference [Pa]')
plt.title('Algal-toxin minitank, pressure difference [Pa]')
plt.savefig('minitank_dp.png', dpi=100)
plt.close(fig)

dp[dp < 0] = 0
u = pitot_velocity(dp, rho=1.15)

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111, xlim=(t0, t1), ylim=(0, 14))
plt.plot(time, u, lw=0.4)
plt.plot(time, running_mean(u, 50), 'r-', lw=1)
plt.grid()
plt.xlabel('Time')
plt.ylabel('Velocity [m/s]')
plt.title('Algal-toxin minitank, velocity [m/s]')
plt.savefig('minitank_u.png', dpi=100)
plt.close(fig)
