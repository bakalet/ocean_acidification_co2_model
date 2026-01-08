# Ocean_acidification_CO2_model

## Authors
- KALETKA Baptiste  
- BARAER Soizic

## Description
A small model to simulate the acidification of a two-layer ocean driven by the increase in atmospheric CO₂ through chemical reactions.

## Databases

### Ocean pH measurements
- `global_omi_health_carbon_ph_area_averaged_1985_P20230930.nc`

### Surface atmospheric CO₂ measurements
- `co2_mlo_surface-insitu_1_ccgg_MonthlyData.nc`

## Code

### Simulation
The two scripts are similar; the difference lies in the parameter optimization method used (`minimize` or `curve_fit`).  
Both are fully functional and independent.

For proper execution, the folder must be extracted and the scripts run from the same directory as the datasets (otherwise, the file paths must be updated in the code).

### Test_opt
Trying to optimize the reaction rate constants (**not functional yet**).

## References

[1] Declining oxygen in the global ocean and coastal waters, Breitburg, D.; Levin, L.A.; et al., *Science*, 359 (6371), 2018, pp. 46  

[2] R.F. Weiss, *Carbon Dioxide in Water and Seawater: The Solubility of a Non-Ideal Gas*, Marine Chemistry, 2 (1974), 203–215  

[3] Keeling, R.F.; Körtzinger, A.; Gruber, N., *Annual Review of Marine Science*, 2 (2010), pp. 199–229  

[4] Li, C.Y.; Huang, J.P.; et al., *Geophysical Research Letters*, 47 (11), 2020  

[5] Gibbs-SeaWater (GSW) Oceanographic Toolbox — TEOS-10  
http://www.teos-10.org  

[6] Encyclopaedia Britannica — State diagram of carbonates  
https://cdn.britannica.com/82/152182-050-719F9121/state-diagram-carbonates-conditions-oceans-2100.jpg  

[7] Mucci, A., *The solubility of calcite and aragonite in seawater*, American Journal of Science, 283 (7), 1983  

[8] Brečević, L.; Nielsen, A.E., *Journal of Crystal Growth*, 98 (3), 1989  

[9] Plummer, L.N.; Busenberg, E., *The solubilities of calcite, aragonite and vaterite in CO₂-H₂O solutions*,  

[10] Morse, J.W.; Millero, F.J., *Geochimica et Cosmochimica Acta*, 44 (1), 1980  

[11] Ryther, J.H., *Photosynthesis in the Ocean as a Function of Light Intensity*, Woods Hole Oceanographic Institution  

[12] Frenette, J.; Demers, S.; Legendre, L.; Dodson, J., *Limnology and Oceanography*, 38 (3), 1993  

[13] Flynn, K.J., *Progress in Oceanography*, 56 (2003), 249–279  

[14] Omrani, S. et al., *Journal of Molecular Liquids*, 345 (2022)  

[15] Li, Y-H., *Geochimica et Cosmochimica Acta*, 38 (1974), 703–714  

[16] Lueker, T.J.; Dickson, A.G.; Keeling, C.D., *Marine Chemistry*, 70 (2000), 105–119  

[17] Sulpis, O. et al., *Nature Geoscience*, 14 (2021), 423–428  

[18] Krumgalz, B.S., *Oceanologica Acta*, 5 (1), 121–128
