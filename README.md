# Hazelnut-project

These files are associated with the Frontiers in Environmental Archaeology article: "Carbon isotope values of hazelnut shells: a new proxy for canopy density".

Authors: Amy K. Styring, Elin Jirdén, Per Lagerås, Mikael Larsson, Arne Sjöström, Karl Ljung

This repository contains:
1. R script for calculating accuracy and precision of measured isotope values in Supplementary Table 3.
2. Supplementary Table 3: Raw and calibrated isotope data for modern and archaeological hazelnut samples and reference materials.
3. R script for statistical analysis in Styring et al. Carbon isotope values of hazelnut shells: a new proxy for canopy density. Frontiers in Environmental Archaeology.
4. Styring et al_HazelnutData: Data to be used with R script for statistical analysis.

Below are explanatory notes for column headings in Styring_et_al_HazelnutData.xlsx

1. SampleID: Unique Sample ID.
2. Site: Modern sampling site or archaeological site.
3. Location: Location (sub-site) within modern sampling site where a leaf area index measurement was taken and hazelnuts sampled from individual trees or from the ground below trees.
4. Tree: Where nuts were collected directly from a tree, or there was only one hazel tree in the location, a unique Tree ID is given.
5. Intra-tree/nut: Indicates which modern nutshell samples were included in the intra-tree or intra-nutshell comparisons.
6. LAI: Leaf Area Index measured at each modern sampling location.
7. d13CCO2: The d13C value of atmospheric carbon dioxide. For modern samples this is an annual average d13C value determined from air sampled at weekly intervals during 2021 and 2021 (from Global Monitoring Laboratory https://gml.noaa.gov/dv/iadv/, Pallas-Sammaltunturi station, Finland). For archaeological samples this is the d13C value approximated by the AIRCO2_LOESS system (Ferrio et al., 2005).
8. Period: Broad time period from which samples derive.
9. Date: Date or date range of each sample.
10. Acid-treated: Yes (Y) or no (N). See section 2.3 in paper for protocol.
11. Charred: Yes (Y) or no (N).
12. NugR: Raw value of micrograms (ug) of nitrogen measured.
13. d15NR: Raw d15N value measured.
14. CugR: Raw value of micrograms (ug) of nitrogen measured.
15. d13CR: Raw d13C value measured.
16. d18OR: Raw d18O value measured (not used).
17. Nugdc: Micrograms (ug) of nitrogen measured after drift correction.
18. d15Ndc: d15N value after drift correction.
19. Cugdc: Micrograms (ug) of carbon measured after drift correction.
20. d13Cdc: d13C value after drift correction.
21. d18Odc: d18O value after drift correction (not used).
22. Runfile: Date of isotope run.
23. pcC: %C
24. pcN: %N
25. CN: C:N atomic ratio
26. normd13C: d13C value after drift correction, normalised (calibrated) to the VPDB scale using two reference materials.
27. normd15N: d15N value after drift correction, normalised (calibrated) to the AIR scale using two reference materials.
28. D13C: D13C value, calculated from the normdd13C and d13CCO2 values using the equation in section 2.4 of the paper.
29. LAI_bin: The category of canopy density assigned based on measured LAI values. Details are in section 2.4 of the paper.
