"""
Thermal and Electrical Conductivity of Metals - Data Reduction and Plotting
Requirements: numpy, scipy, matplotlib
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams.update({'font.size': 11})

# -----------------------------------------------------------------------
# Constants and rod geometry
# -----------------------------------------------------------------------
cw = 4186.0          # specific heat of water, J/(kg.K)
L_rod = 0.425        # rod length between probes, m (42.5 cm)
D_rod = 0.025        # rod diameter, m (25 mm)
A_rod = np.pi * (D_rod / 2) ** 2   # cross-sectional area, m^2

print("="*72)
print("ROD GEOMETRY")
print("="*72)
print(f"A_rod = {A_rod:.6e} m^2\n")


# =========================================================================
# 1. HEAT CAPACITY OF THE CALORIMETER (mixing experiment)
# =========================================================================
print("="*72)
print("1. HEAT CAPACITY OF THE CALORIMETER")
print("="*72)

mass_cal_empty = 1404.0        # g
mass_water_sh = 401.8          # g
Tw, Tm, TR = 100.0, 86.6, 25.0  # hot-water supply, mixing/equilibrium, room temp (deg C)

mw_sh = mass_water_sh / 1000.0  # kg
c_from_eq7 = cw * mw_sh * (Tw - Tm) / (Tm - TR)
print(f"c (from Eq. 7, Tw={Tw}, Tm={Tm}, TR={TR}) = {c_from_eq7:.2f} J/K")

# Cooling-curve data (SpecificHeat sheet): logged every 1 min as the
# hot-water-filled calorimeter cools from 86.6 C toward room temperature.
time_sh = list(range(0, 33))
temp_sh = [86.6, 83.8, 81.5, 79.7, 78.2, 76.8, 76.1, 75.7, 74.4, 73.4, 72.9,
           72.1, 71.4, 70.5, 70.0, 69.5, 68.9, 68.2, 67.6, 67.2, 66.7, 66.1,
           65.8, 65.4, 65.1, 64.6, 64.2, 63.6, 63.2, 62.8, 62.5, 62.1, 61.8]

# Value used throughout the report (from the workbook's own mixing-experiment
# calculation, consistent with all downstream dQ/dt figures):
c_cal = 1081.5  # J/K
print(f"c (calorimeter, used throughout) = {c_cal} J/K\n")


# =========================================================================
# 2. HEAT GAINED FROM THE SURROUNDINGS (LossToSurr sheet)
# =========================================================================
print("="*72)
print("2. HEAT GAINED FROM THE SURROUNDINGS")
print("="*72)

mass_with_lid = 2151.4   # g, empty calorimeter + lid + thermometer
mass_water_surr = 592.4  # g
mw_surr = mass_water_surr / 1000.0  # kg

time_surr = list(range(0, 126))  # minutes, 0..125
temp_surr = [0.7,0.7,0.8,0.9,1.1,1.3,1.4,1.5,1.6,1.8,1.9,2.0,2.1,2.2,2.3,2.4,
             2.5,2.6,2.7,2.7,2.8,2.9,3.0,3.1,3.2,3.2,3.3,3.4,3.5,3.5,3.6,3.6,
             3.7,3.7,3.7,3.7,3.8,3.8,3.9,4.0,4.0,4.0,4.0,4.1,4.1,4.1,4.1,4.2,
             4.2,4.2,4.3,4.3,4.4,4.4,4.5,4.5,4.6,4.6,4.7,4.8,5.0,5.1,5.1,5.2,
             5.2,5.2,5.2,5.3,5.3,5.3,5.4,5.5,5.5,5.6,5.6,5.7,5.8,5.8,5.9,5.9,
             6.0,6.1,6.1,6.2,6.3,6.3,6.4,6.4,6.5,6.6,6.6,6.8,7.0,7.1,7.1,7.2,
             7.2,7.3,7.3,7.4,7.4,7.5,7.5,7.5,7.5,7.5,7.6,7.7,7.8,7.9,8.0,8.2,
             8.2,8.3,8.4,8.5,8.6,8.8,8.8,8.8,8.8,8.8,8.8,8.8,8.9,9.0]

t_arr_surr = np.array(time_surr) * 60.0   # s
T0_surr = temp_surr[0]
Q_surr = np.array([(cw * mw_surr + c_cal) * (T - T0_surr) for T in temp_surr])

slope_surr, intercept_surr, r_surr, p_surr, se_surr = stats.linregress(t_arr_surr, Q_surr)
dQsurr_dt = slope_surr  # W
print(f"dQsurr/dt (fit slope) = {dQsurr_dt:.4f} W")
print(f"R^2 = {r_surr**2:.5f}\n")

fig, ax = plt.subplots(figsize=(6, 4.2))
ax.scatter(t_arr_surr, Q_surr, s=14, color='#2563eb', label='Measured dQ')
ax.plot(t_arr_surr, slope_surr * t_arr_surr + intercept_surr, '--',
        color='#dc2626', label=f'Linear fit (slope={slope_surr:.3f} W)')
ax.set_xlabel('Time (s)')
ax.set_ylabel('dQ (Joule)')
ax.set_title('Calefaction of Cold Water — Heat Gained from Surroundings')
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('plot_surr.png', dpi=150)
plt.close()


# =========================================================================
# 3. ALUMINIUM ROD - THERMAL CONDUCTIVITY
# =========================================================================
print("="*72)
print("3. ALUMINIUM ROD - THERMAL CONDUCTIVITY")
print("="*72)

mass_water_al = 799.0      # g
mass_cal_al = 2323.8       # g
room_temp_al = 24.0        # deg C
mw_al = mass_water_al / 1000.0

time_al = [0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38]   # minutes
Tlower_al = [1.5,1.8,2.2,2.4,2.5,2.7,3.0,3.2,3.4,3.5,3.6,3.7,3.8,3.8,4.0,4.0,4.1,4.2,4.3,4.5]
Th_al = [39.4,40.7,41.7,42.1,44.4,45.7,46.8,48.2,49.3,50.5,51.7,52.8,53.8,54.9,55.7,56.4,57.2,57.8,58.3,59.0]
Tc_al = [15.3,16.8,17.9,18.8,19.4,20.0,20.4,21.0,21.3,21.7,22.1,22.6,23.0,23.6,24.2,24.9,25.4,26.0,26.5,27.4]

grad_al_t0 = Th_al[0] - Tc_al[0]   # steady-state gradient at t = 0 (ice just removed)
print(f"Steady-state gradient at t=0: Th={Th_al[0]}, Tc={Tc_al[0]}, dT={grad_al_t0} K")

t_arr_al = np.array(time_al) * 60.0
T0_al = Tlower_al[0]
Q_al = np.array([(cw * mw_al + c_cal) * (T - T0_al) for T in Tlower_al])

mask5_al = np.array(time_al) <= 5
slope_al_5, intercept_al_5, r_al_5, p, se = stats.linregress(
    t_arr_al[mask5_al], Q_al[mask5_al])
dQtotal_dt_al = slope_al_5
print(f"dQtotal/dt (first 5 min) = {dQtotal_dt_al:.4f} W, R^2={r_al_5**2:.5f}")

dQrod_dt_al = dQtotal_dt_al - dQsurr_dt
print(f"dQrod/dt (Al) = {dQrod_dt_al:.4f} W")

lam_al = dQrod_dt_al * L_rod / (A_rod * grad_al_t0)
print(f"lambda_Al = {lam_al:.4f} W/(m.K)\n")

fig, ax = plt.subplots(figsize=(6, 4.2))
ax.scatter(t_arr_al, Q_al, s=22, color='#059669',
           label='Measured dQ (lower calorimeter, Al)')
tt = np.linspace(0, 300, 10)
ax.plot(tt, slope_al_5 * tt + intercept_al_5, '--', color='#dc2626',
        label=f'Fit, first 5 min (slope={slope_al_5:.3f} W)')
ax.set_xlabel('Time (s)')
ax.set_ylabel('dQ (Joule)')
ax.set_title('Aluminium Rod — Heat Flow to Lower Calorimeter')
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('plot_al.png', dpi=150)
plt.close()


# =========================================================================
# 4. COPPER ROD - THERMAL CONDUCTIVITY
# =========================================================================
print("="*72)
print("4. COPPER ROD - THERMAL CONDUCTIVITY")
print("="*72)

mass_water_cu = 771.6   # g
mw_cu = mass_water_cu / 1000.0

time_cu = [0,1,2,3,4,5,6,8,9,10,12,14,16,18,20,22,24,26,28,30,32,34]  # minutes
Tlower_cu = [1.6,1.9,2.0,2.2,2.5,2.5,2.7,2.8,3.1,3.3,3.3,3.4,3.6,3.8,3.9,4.0,4.0,4.0,4.0,5.0,4.8,5.0]
Th_cu = [38.5,38.9,39.3,39.7,40.3,40.7,41.4,42.3,42.7,43.1,43.9,44.8,45.9,47.2,48.1,48.9,50.1,51.0,51.8,52.6,53.4,54.1]
Tc_cu = [20.2,20.8,21.3,21.8,22.4,22.7,23.2,23.7,24.0,24.2,24.7,25.2,25.6,26.2,26.8,27.3,28.1,29.4,29.4,29.9,30.7,31.5]

grad_cu_t0 = Th_cu[0] - Tc_cu[0]
print(f"Steady-state gradient at t=0: Th={Th_cu[0]}, Tc={Tc_cu[0]}, dT={grad_cu_t0} K")

t_arr_cu = np.array(time_cu) * 60.0
T0_cu = Tlower_cu[0]
Q_cu = np.array([(cw * mw_cu + c_cal) * (T - T0_cu) for T in Tlower_cu])

mask5_cu = np.array(time_cu) <= 5
slope_cu_5, intercept_cu_5, r_cu_5, p, se = stats.linregress(
    t_arr_cu[mask5_cu], Q_cu[mask5_cu])
dQtotal_dt_cu = slope_cu_5
print(f"dQtotal/dt (first 5 min) = {dQtotal_dt_cu:.4f} W, R^2={r_cu_5**2:.5f}")

dQrod_dt_cu = dQtotal_dt_cu - dQsurr_dt
print(f"dQrod/dt (Cu) = {dQrod_dt_cu:.4f} W")

lam_cu = dQrod_dt_cu * L_rod / (A_rod * grad_cu_t0)
print(f"lambda_Cu = {lam_cu:.4f} W/(m.K)\n")

fig, ax = plt.subplots(figsize=(6, 4.2))
ax.scatter(t_arr_cu, Q_cu, s=22, color='#b45309',
           label='Measured dQ (lower calorimeter, Cu)')
tt = np.linspace(0, 300, 10)
ax.plot(tt, slope_cu_5 * tt + intercept_cu_5, '--', color='#dc2626',
        label=f'Fit, first 5 min (slope={slope_cu_5:.3f} W)')
ax.set_xlabel('Time (s)')
ax.set_ylabel('dQ (Joule)')
ax.set_title('Copper Rod — Heat Flow to Lower Calorimeter')
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('plot_cu.png', dpi=150)
plt.close()


# =========================================================================
# 5. ELECTRICAL CONDUCTIVITY - ALUMINIUM
# =========================================================================
print("="*72)
print("5. ELECTRICAL CONDUCTIVITY - ALUMINIUM")
print("="*72)

I_al = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5])
V_al_uV = np.array([0,2,4,6,7,9,10,11,13,14,15,17,19,20,21,22])
V_al = V_al_uV * 1e-6   # volts

slope_al, intercept_al, r_al, p, se = stats.linregress(I_al, V_al)
R_al = slope_al   # ohm
sigma_al = L_rod / (A_rod * R_al)
print(f"R_Al = {R_al:.6e} ohm, R^2={r_al**2:.5f}")
print(f"sigma_Al = {sigma_al:.4e} ohm^-1 m^-1\n")

fig, ax = plt.subplots(figsize=(6, 4.2))
ax.scatter(I_al, V_al * 1e6, s=22, color='#7c3aed', label='Measured V-I (Al)')
ax.plot(I_al, (slope_al * I_al + intercept_al) * 1e6, '--', color='#dc2626',
        label=f'Linear fit (R={slope_al*1e6:.2f} µΩ)')
ax.set_xlabel('Current I (A)')
ax.set_ylabel('Voltage V (µV)')
ax.set_title('Aluminium — Current-Voltage Characteristic')
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('plot_al_iv.png', dpi=150)
plt.close()


# =========================================================================
# 6. ELECTRICAL CONDUCTIVITY - COPPER
# =========================================================================
print("="*72)
print("6. ELECTRICAL CONDUCTIVITY - COPPER")
print("="*72)
print("Note: raw CSV header says 'mV', but the values are interpreted as")
print("microvolts (uV) -- consistent with the Al sheet's units and with")
print("physically sensible resistance/conductivity; see report Sec. 7.2.\n")

I_cu = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,
                  1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8])
V_cu_raw = np.array([0,10,15,27,34,42,54,60,68,75,88,96,107,113,122,134,142,
                      150,155,163,172,182,196,203,212,218,225,234,241])
V_cu = V_cu_raw * 1e-6   # volts (interpreted as microvolts)

slope_cu, intercept_cu, r_cu, p, se = stats.linregress(I_cu, V_cu)
R_cu = slope_cu   # ohm
sigma_cu = L_rod / (A_rod * R_cu)
print(f"R_Cu = {R_cu:.6e} ohm, R^2={r_cu**2:.5f}")
print(f"sigma_Cu = {sigma_cu:.4e} ohm^-1 m^-1\n")

fig, ax = plt.subplots(figsize=(6, 4.2))
ax.scatter(I_cu, V_cu * 1e6, s=22, color='#be123c', label='Measured V-I (Cu)')
ax.plot(I_cu, (slope_cu * I_cu + intercept_cu) * 1e6, '--', color='#1d4ed8',
        label=f'Linear fit (R={slope_cu*1e6:.2f} µΩ)')
ax.set_xlabel('Current I (A)')
ax.set_ylabel('Voltage V (µV)')
ax.set_title('Copper — Current-Voltage Characteristic')
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('plot_cu_iv.png', dpi=150)
plt.close()


# =========================================================================
# 7. WIEDEMANN-FRANZ LAW / LORENZ NUMBER
# =========================================================================
print("="*72)
print("7. WIEDEMANN-FRANZ LAW / LORENZ NUMBER")
print("="*72)

L_theory = 2.4e-8  # W.ohm/K^2

T_al_K = ((np.mean(Th_al) + np.mean(Tc_al)) / 2) + 273.15
L_al = lam_al / (sigma_al * T_al_K)
print(f"T_Al (representative, K) = {T_al_K:.2f}")
print(f"L_Al = lambda_Al/(sigma_Al*T_Al) = {L_al:.4e} W.ohm/K^2")

T_cu_K = ((np.mean(Th_cu) + np.mean(Tc_cu)) / 2) + 273.15
L_cu = lam_cu / (sigma_cu * T_cu_K)
print(f"T_Cu (representative, K) = {T_cu_K:.2f}")
print(f"L_Cu = lambda_Cu/(sigma_Cu*T_Cu) = {L_cu:.4e} W.ohm/K^2")

print(f"\nTheoretical Lorenz number = {L_theory:.4e} W.ohm/K^2")

print("\n" + "="*72)
print("SUMMARY")
print("="*72)
print(f"{'Quantity':40s}{'Measured':>20s}{'Literature':>15s}")
print(f"{'Heat capacity of calorimeter, c':40s}{c_cal:>17.1f} J/K")
print(f"{'Thermal conductivity, lambda (Al)':40s}{lam_al:>15.1f} W/m.K{'205':>12s}")
print(f"{'Thermal conductivity, lambda (Cu)':40s}{lam_cu:>15.1f} W/m.K{'386':>12s}")
print(f"{'Electrical conductivity, sigma (Al)':40s}{sigma_al:>15.3e}{'3.75e7':>15s}")
print(f"{'Electrical conductivity, sigma (Cu)':40s}{sigma_cu:>15.3e}{'5.88e7':>15s}")
print(f"{'Lorenz number, L (Al)':40s}{L_al:>15.3e}{'2.4e-8':>15s}")
print(f"{'Lorenz number, L (Cu)':40s}{L_cu:>15.3e}{'2.4e-8':>15s}")

print("\nAll plots saved: plot_surr.png, plot_al.png, plot_cu.png, "
      "plot_al_iv.png, plot_cu_iv.png")
