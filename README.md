# Sea-level Rise Model
This is the sea-level rise component as it is/will be implemented into the new integrated assessment model [FRIDA](https://github.com/metno/WorldTransFRIDA). It takes global mean surface temperature anomalies since 1850 and ocean heat content changes as input to calculate total global sea-level rise (SLR) as the sum of five contributions:
thermal expansion of the ocean, mountain glacier melt, Greenland ice sheet melt, Antarctic ice sheet melt and land water storage changes.

## Thermosteric SLR
The annual controibution of thermosteric SLR is calculated as a linear function of ocean heat content (OHC) changes, with a default factor of 0.105 (m / YJ).
```math
\frac{dSLR_{thermo}}{dt} = r \cdot \frac{dOHC}{dt}
```
Several different values for $r$ can be found in the literature:
- 0.145 (Marti et al., 2023)
- 0.12 (Levitus et al., 2012)
- 0.15 (Church et al., 2011)
- 0.105 (BRICKv0.3; Wong et al., 2017)
- 0.115 from fitting MPI-ESM data
- 0.1199 from IPCC estimates of ocean heat content changes and thermosteric SLR
- 0.0975 from fitting observational data of ocean heat content and thermosteric SLR

A suggested uncertainty range here is [0.09, 0.12], as with higher values the fit to observational data (ESA SLBC CCI) is bad, when using OHC input from FaIR.

## SLR from mountain glacier melt
This is calculated using a widely used parameterization (e.g. in BRICKv0.3; Wong et al., 2017) from Wigley and Raper (2005):
```math
\frac{dSLR_{MG}}{dt} = \beta_0\cdot T(t)^{p}\cdot(1.0 - \frac{SLR_{MG}}{V_0})^n
```
The nonlinearity exponent $p$ was added to the temperature anomaly, as this improves the fit with IPCC estimates (Fox-Kemper et al., 2021: Table 9.9). The parameter values are chosen as in Perette et al. (2013) and Li et al. (2020).
$p$ is set to 1.5. The uncertainty within the input global mean temperature leads to enough uncertainty of future SLR from this component compared to IPCC estimates (a consequence of $p>1$) that no additional range of parameter uncertainty needs to be given.

## SLR from Greenland ice sheet
The SLR contribution of the Greenland ice sheet is calculated as in MAGICC6 (Nauels et al., 2017), dividing the contribution into changes in the surface mass balance and the discharge. A scale factor was added, which is multiplied to three of the parameters within the calculation to allow for a simple expression of uncertainties. $SLR_{GIS}$ is higher, if this scale factor is higher, and values should be within the range [0.7,1.0] (default: 0.85). 

## SLR from Antarctic ice sheet
This implementation follows that of the BRICKv0.3 model (Wong et al., 2017), which employs the Danish Center for Earth System Science Antarctic Ice Sheet (DAIS) model to simulate the Antarctic Ice Sheet contribution to global sea level (Shaffer, 2014).
An option to add Marine Ice Cliff Instability (MICI) behaviour was included in the model, as it is done in [ḾimiBRICK.jl](https://github.com/raddleverse/MimiBRICK.jl?tab=readme-ov-file).
In order to add some expression of uncertainty to this component, a parameter called "AIS_anto_addend" with possible values between 0.02 and 0.07 (default 0.04) was introduced, which is added to the two parameters used in the calculation of Antarctic subsurface ocean temperature.

## SLR from Land water storage
In FRIDA this currently is a function of population. Here, it is simply a constant rate (0.3 mm/year) plus some noise.

## References:
- Church, J. A., White, N. J., Konikow, L. F., Domingues, C. M., Cogley, J. G., Rignot, E., Gregory, J. M., van den Broeke, M. R., Monaghan, A. J., and Velicogna, I.: Revisiting the Earth's sea-level and energy budgets from 1961 to 2008, Geophys. Res. Lett., 38, L18601, https://doi.org/10.1029/2011GL048794, 2011.
- Fox-Kemper, B., H.T. Hewitt, C. Xiao, G. Aðalgeirsdóttir, S.S. Drijfhout, T.L. Edwards, N.R. Golledge, M. Hemer, R.E. Kopp, G.  Krinner, A. Mix, D. Notz, S. Nowicki, I.S. Nurhati, L. Ruiz, J.-B. Sallée, A.B.A. Slangen, and Y. Yu, 2021: Ocean, Cryosphere and Sea Level Change. In Climate Change 2021: The Physical Science Basis. Contribution of Working Group I to the Sixth Assessment Report of the Intergovernmental Panel on Climate Change: Masson-Delmotte, V., P. Zhai, A. Pirani, S.L.  Connors, C. Péan, S. Berger, N. Caud, Y. Chen, L. Goldfarb, M.I. Gomis, M. Huang, K. Leitzell, E. Lonnoy, J.B.R. Matthews, T.K. Maycock, T. Waterfield, O. Yelekçi, R. Yu, and B. Zhou (eds.). Cambridge University Press, Cambridge, United Kingdom and New York, NY, USA, pp. 1211–1362, https://doi.org/10.1017/9781009157896.011.
- Levitus, S., Antonov, J. I., Boyer, T. P., Baranova, O. K., Garcia, H. E., Locarnini, R. A., Mishonov, A. V., Reagan, J. R., Seidov, D., Yarosh, E. S., and Zweng, M. M.: World ocean heat content and thermosteric sea level change (0–2000 m), 1955–2010, Geophys. Res. Lett., 39, L10603, https://doi.org/10.1029/2012GL051106, 2012. 
- Li et al., Optimal temperature overshoot profile found by limiting global sea level rise as a lower-cost climate target.Sci. Adv.6,eaaw9490(2020). https://doi.org/10.1126/sciadv.aaw9490
- Marti, F., Blazquez, A., Meyssignac, B., Ablain, M., Barnoud, A., Fraudeau, R., Jugier, R., Chenal, J., Larnicol, G., Pfeffer, J., Restano, M., and Benveniste, J.: Monitoring the ocean heat content change and the Earth energy imbalance from space altimetry and space gravimetry, Earth Syst. Sci. Data, 14, 229–249, https://doi.org/10.5194/essd-14-229-2022, 2022.
- Nauels, A., Meinshausen, M., Mengel, M., Lorbacher, K., and Wigley, T. M. L.: Synthesizing long-term sea level rise projections – the MAGICC sea level model v2.0, Geosci. Model Dev., 10, 2495–2524, https://doi.org/10.5194/gmd-10-2495-2017, 2017.
- Perrette, M., Landerer, F., Riva, R., Frieler, K., and Meinshausen, M.: A scaling approach to project regional sea level rise and its uncertainties, Earth Syst. Dynam., 4, 11–29, https://doi.org/10.5194/esd-4-11-2013, 2013.
- Shaffer, G.: Formulation, calibration and validation of the DAIS model (version 1), a simple Antarctic ice sheet model sensitive to variations of sea level and ocean subsurface temperature, Geosci.Model Dev., 7, 1803–1818, https://doi.org/10.5194/gmd-7-1803-2014, 2014.
- Wigley, T. M. L., and S. C. B. Raper (2005), Extended scenarios for glacier melt due to anthropogenic forcing, Geophys. Res. Lett., 32, L05704, https://doi.org/10.1029/2004GL021238.
- Wong, T. E., Bakker, A. M. R., Ruckert, K., Applegate, P., Slangen, A. B. A., and Keller, K.: BRICK v0.2, a simple, accessible, and transparent model framework for climate and regional sea-level projections, Geosci. Model Dev., 10, 2741–2760, https://doi.org/10.5194/gmd-10-2741-2017, 2017
