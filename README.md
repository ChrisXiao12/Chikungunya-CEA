The code is provided for free use in research and educational purposes.
The authors assume no responsibility or liability for any errors, omissions, or outcomes resulting from its use. Users assume full responsibility for verifying results and ensuring suitability for their specific applications.

Suggested citation: Drabo, E., & Xiao, C. (2026). Chikungunya Vaccine SVEIRD Model (R implementation). Available at: https://github.com/ChrisXiao12/Chikungunya-CEA

Instructions
Enter CHIKV_CEA.R and set directories as appropriate. Set simulation settings as desired.

Intermediate RDS outputs are stored as GitHub LFS due to size. 

Flow_diagram creates the SVEIRD model. Functions loads the needed functions. Data inputs the data from parameters.xlsx.
Analyses runs the base-case, one-way sensitivity analysis, additional sensitivity analyses, and probabilistic sensitivity analyses.

Run-time estimate (16gb ram, AMD Ryzen 9 8945HS processor)
Base-Case: 1 minute
OWSA: 14 minutes
Additional Sensitivity Analyses: ~ 3 hours 40 minutes
PSA: ~ 3 hours and 15 minutes (1000 iterations)
