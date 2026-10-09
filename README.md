# PCB cell culture system modeling: Equilibrium approach

----------------------
General Information
----------------------

Deposit Title: PCB partitioning in multiwell setup experiments

Contributor information:

Andres Martinez, PhD
University of Iowa - Department of Civil & Environmental Engineering
Iowa Superfund Research Program (ISRP)
andres-martinez@uiowa.edu
ORCID: 0000-0002-0572-1494

This README file was generated on Otober 9 2026 by Andres Martinez.

This work was supported by the National Institutes of Environmental Health Sciences (NIEHS) grant #P42ES013661.  The funding sponsor did not have any role in study design; in collection, analysis, and/or interpretation of data; and/or in the decision to submit this data for publication or deposit it in a repository.

This README file describes a revised version of scripts for calculating the individual fraction of PCBs and their metabolites in multiwell experiments and the creation of Figure 7 in:

"Structure Activity Relationship of Lower Chlorinated Biphenyls and their Human-relevant Metabolites for Astrocyte Toxicity" 

Neha Paranjape1,2, Laura E. Dean2, Andres Martinez3, Hans-Joachim Lehmler2, Jonathan A. Doorn1, * 

1Department of Pharmaceutical Sciences & Experimental Therapeutics, College of Pharmacy, University of Iowa, Iowa City, IA, USA
2Department of Occupational and Environmental Health, College of Public Health, University of Iowa, Iowa City, IA, USA
3Department of Civil and Environmental Engineering, IIHR-Hydroscience & Engineering, University of Iowa, Iowa City, IA, USA
4Department of Environmental and Radiological Health Sciences, College of Veterinary Medicine and Biomedical Sciences, Colorado State University, Fort Collins, CO USA. 

*Corresponding Author:  
"Jonathan A. Doorn 
Professor and Chair 
Department of Pharmaceutical Sciences & Experimental Therapeutics 
College of Pharmacy, 536 CPB 
The University of Iowa 
180 S. Grand Ave. 
Iowa City, IA  52242 USA""
 
Email:jonathan-doorn@uiowa.edu  
TEL: 319-335-8834

The script is developed for individual PCBs, OH-PCBs and sulfate-PCBs from experiments from Project 1, ISRP.

--------
PREREQUISITES & DEPENDENCIES
--------

This section of the ReadMe file lists the necessary software required to run codes in "R".

Software:
- Any web browser (e.g., Google Chrome, Microsoft Edge, Mozilla Firefox, etc.)
- R-studio for easily viewing, editing, and executing "R" code as a regular "R script" file:
https://www.rstudio.com/products/rstudio/download/

--------
SOFTWARE INSTALLATION
--------

This section of the ReadMe file provides short instructions on how to download and install "R Studio".  "R Studio" is an open source (no product license required) integrated development environment (IDE) for "R" and completely free to use.  To install "R Studio" follow the instructions below:

1. Visit the following web address: https://www.rstudio.com/products/rstudio/download/
2. Click the "download" button beneath RStudio Desktop
3. Click the button beneath "Download RStudio Desktop".  This will download the correct installation file based on the operating system detected.
4. Run the installation file and follow on-screen instructions. 

--------
HOW TO RUN THE CODE
--------

It is recommended to create a project in RStudio and name it Metab-Partitioning. Paste and execute Subfolders.R code to generate the different folders. Then paste the MetabPartitioning.R into the R folder just generated. All the packages needed, with their respecting libraries, are included in the code. All the needed data are in the code. If needed, the fraction results can be exported as a csv file to the Output/Data/csv created folder.

--------
EXTERNAL DATA
--------

Partitioning coefficients were obtained from UFZ-LSER database:

Ulrich, N., Endo, S., Brown, T.N., Watanabe, N., Bronner, G., Abraham, M.H., Goss, K.-U., UFZ-LSER database v 3.2.1 [Internet], Leipzig, Germany, Helmholtz Centre for Environmental Research-UFZ. 2017 [accessed on 08.05.2023]. Available from http://www.ufz.de/lserd.
Metabolites were obtained using their SMILE configuration.
