# OmniPath-based generation of the NicheNet Method prior model: A SARS-CoV-2 case study

## Overview

NicheNet is a recently developed method to prioritize ligand–target 
relationships between interacting cells by combining their expression data 
with prior knowledge on interaction networks. For this purpose, it explores the 
most consistent inter- and intra-cellular protein interactions in accordance 
with a given gene expression dataset. The authors collected different types of 
interactions from more than 20 databases to build a ligand-receptor network, a 
signaling network and a gene regulatory network. OmniPath provides a 
single-access point covering all the different types of interactions employed in 
the NicheNet method. **Therefore, we here highlight the value of OmniPath by 
exclusively using its resources to create a ligand-target regulatory potential 
model as described in the NicheNet article** (Browaeys, Saelens and Saeys, 2019).     

Nowadays, COVID-19, caused by SARS-CoV-2, is spreading globally throughout the 
planet. WHO has reported  approximately 10 million confirmed cases and 500 000 
deaths to date (June 29, 2020). Against this background, we aim at using our 
Omnipath-based version of NicheNet to explore the autocrine signaling after 
SARS-CoV-2 infection. In particular, **we explore the potential regulatory 
effect of over-expressed ligands after infection on the expression of 
inflammatory response related genes in the Calu3 cell line**. The RNAseq 
expression data was taken from a recent study. 

## Content

This repository contains the vignettes to reproduce the SARS-CoV-2 case study
presented in the publication: 

> Türei, D., Valdeolivas, A., Gul, L.,  Palacio-Escat, N., Ivanova, O., Modos, D., Korcsmáros T. & Saez-Rodriguez, J. 
Integration of intra- and intercellular signaling resources with OmniPath



## Important Links

### Omnipath:
<http://omnipathdb.org/>
<br>
<https://github.com/saezlab/pypath>
<br>
<https://saezlab.github.io/OmnipathR/>

### NicheNet:
<https://www.nature.com/articles/s41592-019-0667-5>
<br>
<https://github.com/saeyslab/nichenetr/>

### RNAseq Expression data 
<https://www.biorxiv.org/content/10.1101/2020.03.24.004655v1>
<br>
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507>

## References

> Türei, D., Korcsmáros, T., & Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. _Nature methods_, 13(12), 966–967. [10.1038/nmeth.4077](https://doi.org/10.1038/nmeth.4077)

> Browaeys, R., Saelens, W. & Saeys, Y. NicheNet: modeling intercellular communication by linking ligands to target genes. _Nature Methods_ 17, 159–162 (2020). [10.1038/s41592-019-0667-5](https://doi.org/10.1038/s41592-019-0667-5)

> Blanco-Melo, D., Nilsson-Payant, B.E., Liu, W.-C., Uhl, S., Hoagland, D., Møller, R., Jordan, T.X., Oishi, K., Panis, M., Sachs, D., Wang, T.T., Schwartz, R.E., Lim, J.K., Albrecht, R.A., tenOever, B.R., 2020. Imbalanced Host Response to SARS-CoV-2 Drives Development of COVID-19. _Cell_ 181, 1036–1045.e9. [10.1016/j.cell.2020.04.026](https://doi.org/10.1016/j.cell.2020.04.026)
