# SCNA Analysis
Somatic Copy Number Alteration analysis with HMF tools.

The analysis is performed using some Hartwig Medical Foundation (HMF) tools and GRIPSS. All steps are connected and controlled using Nextflow workflow manager. 
Steps:
* Input preparation: screening of the samples already analysed or not and creation of suitable lists for each process.
* GRIDSS: structural variants (SVs) caller which takes BAM files and performs a joint calling extracting breakpoints and single breakends at patient level
* PoN: Panel of Normals (PoN) updating. PoN from [HMF resources](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/dna_pipeline;tab=objects?prefix=&forceOnObjectsSortingFiltering=false) is updated with normals present in the used cohort through a GRIDSS function. 
* Gripss: Somatic filtering of the variants using the PoN produced in the previous step. Done at sample level.
* Amber: extraction of the B Allele Frequency (BAF) given a VCF with heterozygous sites in the genome. Done at sample level.
* Cobalt: calculation of the read depth adn GC normalization step. Done at sample level
* Purple: combines data from Gripss, Amber and Cobalt with also SNVs (from another analysis) and calculates ploidy and purity. It performs the final segmentation inferring copy number regions. Given the Driver genes list ([HMF resources](https://console.cloud.google.com/storage/browser/hmf-public/HMFtools-Resources/dna_pipeline;tab=objects?prefix=&forceOnObjectsSortingFiltering=false)) it also tries to understand the driver mutations. Done at sample level.
* Results: production of sunrises plots, multiploidy plots for multiploidy patients. 

![Pipeline figure](pipeline_schema.pdf)
