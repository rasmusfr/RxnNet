This README provides information about the files in rxnnet/data/database-20.03.24/datasets

nomenclature:
hf_rf - predicted enthalpy of formation from Gong et al.; reference doi: https://10.1021/jacsau.2c00235
hf_mppred - predicted enthalpy of formation from MPpredictor; reference doi: https://doi.org/10.1021/acs.jcim.3c00307
e_r2scan - Materials Project r2SCAN electronic energy; reference doi: https://doi.org/10.1103/PhysRevMaterials.6.013801
e_gga_gga_u - Materials Project GGA/GGA+U electronic energy; reference doi: https://doi.org/10.1103/PhysRevB.84.045115
hf_gga_gga_u_r2scan - Materials Project GGA/GGA+U/r2SCAN enthalpy of formation; reference doi: https://doi.org/10.1103/PhysRevMaterials.6.013801
e_chgnet - electronic energy calculated using CHGNet; reference doi: https://doi.org/10.1063/1.555845
e_wang - experimental enthalpy of formation; https://doi.org/10.1038/s41598-021-94550-5
h_ref - NBS tables experimental enthalpy of formation; https://doi.org/10.1063/1.555845


intermetallics.pkl.gz
-pandas dataframe with data pertaining to intermetallic compounds.

lit_disag.pkl.gz
-pandas dataframe with data pertaining to compounds with disagreement between experimental values the NBS data and different sources.

main.pkl.gz
-pandas dataframe with data pertaining to the compounds used in the main text of the paper.

wang_disag.pkl.gz
-pandas dataframe with data pertaining to compounds with disagreement between experimental values the NBS data and the Wang data.
