# MetapopulationBrownRot
Code associated to our publication "A metapopulation framework integrating landscape heterogeneity to model an airborne plant pathogen: The case of brown rot of peach in France" (https://doi.org/10.1016/j.agee.2024.108994)

# A metapopulation framework integrating landscape heterogeneity to model an airborne plant pathogen: the case of brown rot of peach in France

## Description of the data and file structure

This script simulates an epidemiological-metapopulation model describing the multiannual spatiotemporal climate-dependent epidemiological dynamics for brown rot of peaches in southern France. In this model, french peach orchards are subdivided into 755 spatial units (or cells), each one characterized in terms of surface, cultivar and climate, and during ripening period (i.e. when the fruit is susceptible) can develop internal infections (modelled as a susceptible-exposed-infected system) as well as infect other units via airborne epidemic connections (via a stochastic extraction based on spore transport matrices). At the beginning of each year, some units may be already be infected by precedent year mummified fruits (all fruits that have been infected) or because of random inoculum. This model is used to assess epidemic risk indices for each unit, i.e. dangerousness (risk of causing secondary infection in other sites) and vulnerability (risk of becoming infected) across the our study area.

The scripts are:

* *Main.R*, which is supposed to be run, which loads the model parameters, run the metapopulation model, compute and plots dangerousness and vulnerability. This scripts call the following two:
* *Ode\_system.R*, which implements the system of Ordinary Differential Equations describing the in-unit SEIM model.
* *Stoch\_mod.R*, which implements the stochastic part of the model, i.e. external contagion (depending on the rate Rtt) and the interannual persistence of inoculum (depending on the mummified fuit *M* and the probability of random persistence *Ptilde0*)

*Main.R* set the following parameters:

* reps: number of stochastic repetitions of the model
* starting\_years: set of the possible starting years to stochastically estimate the epidemic risk indices
* ly: simulation horizon (number of years)
* Dt: 1 day, integration step of the epidemiological model
* V\_t0: initial susceptible fruit load [fruits/m²]
* E\_t0: initial possible exposed fruit load [fruits/m²]
* Ptilde0: probability of having an inoculum at the beginning of a new ripening season due to external sources

Data are:

* the *domain* shapefile, containing the spatial units of the metapopulation and their cultivated surface
* *thetas.RData*, a data.frame with the set of the model parameters theta\_E, theta\_O and theta\_L identified via the Approximate Bayesian Computation procedure
* Ten files called *List\_params\_i.RData* (i = 1, 2..., 10) containing each one a list of 41 lists of parametres, representing each the epidemiological parameters for a given year in the period 1981-2021 (such as the duration of susceptibility period for each unit). Each file correspond to a randomization of peach cultivars in domain. Each list contains:
    \# eta\, i\.e\. the matrix of the temperature\-dependent spore mortality rate;
    \# lambda\, i\.e\. transmission rate \(beta in eq\. 1 in the article\);
    \# rho\, i\.e\. infection\-related abscission rate \(alpha in eq\. 1 in the article\);
    \# sigma\, i\.e\. matrix of the rain\-dependent infection rate;
    \# g\, d\, b\, i\.e\. parameters of the natural abscission rate;
    \# t\_B\, i\.e\. vector of the blooming dates;
    \# t\_PH\, i\.e\. vector of the pit hardening dates;
    \# t\_0\, i\.e\. earliest pit hardening date among cells;
    \# t\_sim\, i\.e\. vectors with the days of the year corresponding to the integration horizon;
    \# cpd\, i\.e\. vector of the identifiers of cells successfully accomplishing phenological requirements;
    \# GDD\_c\_ec\, i\.e\. vector of the Growing Degree Days of the cultivated variety\. 678 = early\, 1026 = midearly\, 1371 = midlate\, 1772 = late
    \# w\_H\, i\.e\. weight of the fruit at harvest time;
* Ten files called *List\_C\_i.RData* (i = 1, 2..., 10) containing each one a list of 41 lists of matrices, representing the daily connectivity matrices **C** (in the article they are called **W**) among units for each year in the period 1981-2021 during susceptibility period. These matrices considers only rain deposition (i.e., **C<sub>ij(d)</sub>** &gt; 0 only if rain occurred in j at day of the year d). Each file correspond to a randomization of peach cultivars in domain.

## Sharing/Access information

Connectivity matrices have been estimated using HYSPLIT Lagrangian trajectory model

* https://www.ready.noaa.gov/HYSPLIT\_traj.php

Other data:

* SAFRAN grid: https://doi.org/10.57745/1PDFNL
* Peach cultivated surfaces by European regions: https://data.europa.eu/data/datasets/00nlsr3zhd3s6picskoxg?locale=en
* Platform to access the RGA data (e.g., stone fruits by municipality): https://agreste.agriculture.gouv.fr/agreste-web/disaron/Pri2105/detail/
* Brown rot incidence observations in France, used to calibrate the model: https://www.data.gouv.fr/fr/datasets/brown-rot-of-peach-severity-data/
* Weather time series from the SICLIMA portal: https://www6.paca.inrae.fr/agroclim/Demande-des-donnees-de-Meteo-France/SAFRAN-SICLIMA
* Phenological model describing peach ripening depending on temperature: https://www.sciencedirect.com/science/article/pii/S0168192322004798
* Climate-dependent epidemiological model describing in-orchard epidemics: https://www.sciencedirect.com/science/article/pii/S0168192321000101

## Code/Software

The code is written in R. Data are stored within a zipped folder, which needs to be unzipped before running the code.
