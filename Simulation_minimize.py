import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import netCDF4 as nc
from datetime import datetime, timedelta
from scipy.optimize import minimize

# Ouvrir le fichier NetCDF
file_path_CO2 = "co2_mlo_surface-insitu_1_ccgg_MonthlyData.nc"  # Chemin du fichier
dataset = nc.Dataset(file_path_CO2, mode='r')

# Extraction des variables nécessaires
time_decimal = dataset.variables['time_decimal'][:]
co2_concentration = dataset.variables['value'][:]

# Identifier la valeur de remplissage (_FillValue)
fill_value = dataset.variables['value']._FillValue

# Filtrer les valeurs invalides et incohérentes (exclure les valeurs de remplissage et négatives)
valid_mask = (co2_concentration != fill_value) & (co2_concentration > 0)

# Appliquer le masque pour ne garder que les valeurs valides
time_decimal_clean = time_decimal[valid_mask]
co2_concentration_clean = co2_concentration[valid_mask]


# Charger le fichier NetCDF
file_path = "global_omi_health_carbon_ph_area_averaged_1985_P20230930.nc"  # Remplace par le chemin de ton fichier
dataset = nc.Dataset(file_path, mode='r')

# Lire les variables
time = dataset.variables['time'][:]  # Temps en jours depuis 1950
ph = dataset.variables['ph'][:]  # pH

# Convertir le temps en années
base_date = datetime(1950, 1, 1)
years = np.array([base_date + timedelta(days=int(t)) for t in time])
years = np.array([date.year for date in years])  # Extraire seulement l'année

# Épaisseurs des zones
e1 = 150  
e2 = 4000 

# Constantes pour CO2 (en mol/l*atm)
a1_CO2 = -58.0931
a2_CO2 = 90.5069
a3_CO2 = 22.2940
B1_CO2 = 0.027766
B2_CO2 = -0.025888
B3_CO2 = 0.0050578

# Constantes générales
T_ocean = 273.15 + 25 
salinity = 35 
T_ocean2 = 2 + 273.15
Ma = 58.44 #Molar mass of the NaCl


# Constante de Henry
Henry_constant = np.exp(
    a1_CO2 + a2_CO2 * (100 / T_ocean) + a3_CO2 * np.log(T_ocean / 100) + 
    salinity * (B1_CO2 + B2_CO2 * (T_ocean / 100) + B3_CO2 * (T_ocean / 100) ** 2))

def Diffusion_coeff_corre1(T,S):
    C = S*1000/(Ma)
    D = (-18.157948*np.exp(-0.05736*C) + 0.068700205361*T - 0.0003876102395346346*(C*0.820561458)*(T**1.46331515077))*1e-9
    return D

D1 = Diffusion_coeff_corre1(T_ocean,salinity) # coeff diffusion zone 1 m^2/s
D2 = Diffusion_coeff_corre1(T_ocean,salinity) # coeff diffusion zone 2 m^2/s
k_mass1 = D1/e1 # Transfer coefficient zone 1 (m/s)
k_mass2 = D2/e2 # mass transfer coeff between zone 1 and 2 (m/s)
k_transf_CO2 = (k_mass1*k_mass2)/(k_mass1+k_mass2) # mass transfer coeff apparent between the 2 zones (m/s)

nb_year = 37
t_max = nb_year * 3.171e8
t_eval = np.linspace(0, t_max, 40)

# Constantes de réaction
pK1_35_25=6.000
pK1_35_2=6.017
pK3_35_25=9.115
pK3_35_2=9.431

k1 = 1.5e-14
k_1 = k1*(10**(-pK1_35_25))
k3 = 1.1e-10
k_3 = k3*(10**(-pK3_35_25))


#Diffusion coefficient m^2/s
D_H_zone1 = 93.1e-10
D_H_zone2 = 56.1e-10

k_H_mass1 = D_H_zone1/e1 # Transfer coefficient zone 1 (m/s)
k_H_mass2 = D_H_zone2/e2 # mass transfer coeff between zone 1 and 2 (m/s)
k_transf_H = (k_H_mass1*k_H_mass2)/(k_H_mass1+k_H_mass2) # mass transfer coeff apparent between the 2 zones (m/s)

D_HCO3_zone1 = 11.8e-10
D_HCO3_zone2 = 2.8e-10 #Estimated 

k_HCO3_mass1 = D_HCO3_zone1/e1 # Transfer coefficient zone 1 (m/s)
k_HCO3_mass2 = D_HCO3_zone2/e2 # mass transfer coeff between zone 1 and 2 (m/s)
k_transf_HCO3 = (k_HCO3_mass1*k_HCO3_mass2)/(k_HCO3_mass1+k_HCO3_mass2) # mass transfer coeff apparent between the 2 zones (m/s)

D_CO32_zone1 = 9.55e-10
D_CO32_zone2 = 4.39e-10 

k_CO32_mass1 = D_CO32_zone1/e1 # Transfer coefficient zone 1 (m/s)
k_CO32_mass2 = D_CO32_zone2/e2 # mass transfer coeff between zone 1 and 2 (m/s)
k_transf_CO32 = (k_CO32_mass1*k_CO32_mass2)/(k_CO32_mass1+k_CO32_mass2) # mass transfer coeff apparent between the 2 zones (m/s)

####### Calcium #########
K_am = -6.1987 - 0.005336 - 0.0001096 * (T_ocean2-273.15)
K_vat = - 172.1295 - 0.077993*T_ocean2 + 3074.688/T_ocean2 + 71.595*np.log(T_ocean2)
K_calc = -171.9065 - 0.077993*T_ocean2 + 2839.319 + 71.595*np.log(T_ocean2)

CO2_ini_ppm=co2_concentration_clean[0]
P_CO2_ini = (1 * CO2_ini_ppm / 1e6)
CO2_ini = P_CO2_ini * Henry_constant 

H_ini = 10**(-ph[0])  # Conversion du pH en [H+]

# Conditions initiales [CO2_zone1, CO2_zone2, C_H_zone1, C_H_zone2, C_HCO3_zone1, C_HCO3_zone2, C_CO32_zone1, C_CO32_zone2]
C_0 = np.array([CO2_ini, 0, H_ini, H_ini, 0, 0, 0, 0, 5.46e-4])

def acc_co2(t, C):
    # Interpolation de la concentration de CO2 atmosphérique
    year_t = 1985 + t / 3.171e8
    co2_atm = np.interp(year_t, time_decimal_clean, co2_concentration_clean)
    P_CO2 = (1 * co2_atm / 1e6)  # Conversion en atm

    # Paramètres
    C_sat = P_CO2 * Henry_constant   # Concentration de saturation (mol/L)

    CO2_zone1, CO2_zone2, C_H_zone1, C_H_zone2, C_HCO3_zone1, C_HCO3_zone2, C_CO32_zone1, C_CO32_zone2, C_Ca2 = C  

    r1 = k1 * CO2_zone1
    r_1 = k_1 * C_H_zone1 * C_HCO3_zone1
    r3 = k3 * C_HCO3_zone1
    r_3 = k_3 * C_CO32_zone1 * C_H_zone1

    dy_dt = np.zeros(9)

    # Calcul de la réaction de dissolution du calcium
    r4 = 6.34e-16  # mol.kg^-1.s^-1
    epsilon = 1e-10  # Pour éviter les divisions par zéro
    
    # Calcul de k4_calc avec une protection contre la division par zéro
    k4_calc = r4 / (max(dy_dt[7] * dy_dt[8], epsilon)) 
    if np.isnan(k4_calc) or np.isinf(k4_calc):
        k4_calc = 0
    k_4_calc = k4_calc / K_calc
    r_4 = k_4_calc * dy_dt[7] * dy_dt[8] if not np.isnan(k_4_calc) and not np.isinf(k_4_calc) else 0

    # CO2 Zone 1
    dy_dt[0] = (-r1 + r_1) + kla * e1 * (C_sat - CO2_zone1) - k_transf_CO2 * (CO2_zone1 - CO2_zone2)
    # CO2 Zone 2
    dy_dt[1] = k_transf_CO2 * (CO2_zone1 - CO2_zone2)
    # H+ zone 1
    dy_dt[2] = (r1 - r_1 + r3 - r_3) - k_transf_H * e1 * (C_H_zone1 - C_H_zone2)
    # H+ zone 2
    dy_dt[3] = k_transf_H * e1 * (C_H_zone1 - C_H_zone2)
    # HCO3- zone 1
    dy_dt[4] = (-r3 + r_3 + r1 - r_1) - k_transf_HCO3 * e1 * (C_HCO3_zone1 - C_HCO3_zone2)
    # HCO3- zone 2
    dy_dt[5] = k_transf_HCO3 * e1 * (C_HCO3_zone1 - C_HCO3_zone2)
    # CO32- zone 1
    dy_dt[6] = (r3 - r_3) - k_transf_CO32 * e1 * (C_CO32_zone1 - C_CO32_zone2)
    # CO32- zone 2
    dy_dt[7] = -r4 + k_transf_CO32 * e1 * (C_CO32_zone1 - C_CO32_zone2)
    #Ca2+ zone 2
    dy_dt[8] = (r_4 - r4)

    return dy_dt

ph_interpolé = np.interp(t_eval / 3.171e8 + 1985, years, ph)

def objective(params):
    global k1, k3, kla
    k1, k3, kla = params
    solution = solve_ivp(acc_co2, [0, t_max], C_0, t_eval=t_eval, method="LSODA")
    pH_simulé = -np.log10(solution.y[2])

    return np.sum((pH_simulé - ph_interpolé) ** 2)

result = minimize(objective, [1e-14, 1e-10, 0.2], method='Nelder-Mead', options={'maxiter': 100})
k1, k3, kla = result.x
print(f"Valeurs optimisées : k1 = {k1:.2e}, k3 = {k3:.2e}, kla = {kla:.2e}")

solution = solve_ivp(acc_co2, [0, t_max], C_0, t_eval=t_eval, method="LSODA")
pH = -np.log10(solution.y[2])

plt.figure(figsize=(10, 5))
plt.plot(years, ph, marker='o', linestyle='-', color='b', label='pH mesuré')
plt.plot(solution.t / 3.171e8 + 1985, pH, label="pH calculé (Zone 1)", color='r', linestyle='-')
plt.xlabel("Année")
plt.ylabel("pH")
plt.title("Évolution du pH de l'eau de mer (Optimisation méthode minimize)")
plt.legend()
plt.grid()
plt.show()