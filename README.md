**MPox-Reckoners: Ready‑Reckoner Modelling Framework for MPox Clade Ib**

A ready‑reckoner modelling tool developed to characterise mpox Clade Ib transmission and control in England. The code accompanies the paper: Targeted vaccination is effective for mpox clade Ib in England despite increased household transmission: predictions from a modelling study by Brooks-Pollock and Danon (2025) https://www.medrxiv.org/content/10.1101/2025.05.27.25328398v1. 

**Summary**
- The ready reckoner framework is applied to Mpox using sexual contacts and household size from National Survey of Sexual Attitudes and Lifestyles 3 (NATSAL3).
- The impact of increasing household secondary attack rates (SAR) was explored, as associated with mpox clade Ib.
- There is a high degree of heterogeneity in individual reproduction numbers. Individuals reporting both same sex and opposite sex sexual contacts disproportionately contributed to transmission potential. 
- There is a minimal impact of increasing household SAR on the overall reproduction number, due to the structure of the sexual contact network and its coincidence with household structure in England.

This is designed as a rapid, network-informed modelling tool to inform strategic outbreak response. This version is calibrated to mpox Clade Ib in England, however the framework is flexible and can be adapted to:
- Other geographical settings informed by local contact survey data.
- Other pathogens with analogous contact types relevant to pathogen-specific transmission routes.
- Updated assumptions on vaccine effectiveness or transmission.

**Repository Contents**

Data sources:
National Survey of Sexual Attitudes and Lifestyles, 2010-2012
NATSAL III; Natsal-3
	
Johnson, A., University College London
Mercer, C., University College London
Sonnenberg, P., University College London
Copas, A., University College London
Wellings, K., London School of Hygiene and Tropical Medicine
Macdowall, W., London School of Hygiene and Tropical Medicine
Erens, B., National Centre for Social Research

https://doi.org/10.5255/UKDA-SN-7799-2


Code: 

**Description and Implementation**
- Individual reproduction numbers (R_i) are calculated using a weighted sum across settings (household, same-sex, opposite-sex contacts), combining contact counts with setting-specific secondary attack rates
- Population-level R_eff is computed by a weighted aggregation of R_i over the population distribution.
- Vaccination reduces R_i according to efficacy against infection and transmission, differing between pre- and post-exposure scenarios.
- The model simulates multiple vaccination scenarios to compare required coverage for control under different strategies and household attack rate assumptions.

To run the model, 

**Output**

**Configuration options**
