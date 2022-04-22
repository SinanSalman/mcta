# Markov Chain Traffic Assignment (MCTA)
This is a Python implementation of the MCTA approach proposed by [1] and expanded upon in [2-5].

The following table describes the use of each mcta module. Examples of its use can be found in the jupyter files `test_mcta.ipynb` and `test_mcta_runtime.ipynb`

File|Description
---|---
mcta_cli.py|MCTA Command Line Interface (CLI). *This is where you should start.*
mcta.py|MCTA module. *This is where the main MCTA logic reside.*
mcta_edit.py|MCTA network modification module
mcta_rw.py|MCTA input/output data read/write module
mcta_vis.py|MCTA visualization module

# Genetic Algorithm
`ga.py` and `ga_batch_run.py` implements a genetic algorithm metaheuristic on top of MCTA. resutls are graphed using `ga_polts.ipynb`.

# emissions.py
A module for estimating emissions based on road network traffic conditions. Estimation methedology is based on: Bureau of Public Roads (BPR) curve, COPERT_v5, and TRANSYT7f.

# Referances
[1] Crisostomi, E., Kirkland, S., & Shorten, R. (2011). A Google-like model of road network dynamics and its application to regulation and control. International Journal of Control, 84(3), 633–651. https://doi.org/10.1080/00207179.2011.568005

[2] Salman, S., & Alaswad, S. (2017). Urban road network crisis response management: Time-sensitive decision optimization. Proceedings of the 2017 Industrial and Systems Engineering Conference, 1307–1313. Retrieved from Scopus. http://amz.xcdsystem.com/iisePapers/iise/2018/papers/SubmitFinalPaper_1350_0304110129.pdf

[3] Salman, S., & Alaswad, S. (2018). Alleviating road network congestion: Traffic pattern optimization using Markov chain traffic assignment. Computers & Operations Research, 99, 191–205. https://doi.org/10.1016/j.cor.2018.06.015

[4] Salman, S., & Alaswad, S. (2019). Mitigating the Impact of Congestion Minimization on Vehicles’ Emissions in a Transportation Network. Proceedings for th 25th International Joint Conference on Industrial Engineering and Operations Management. http://doi.org/10.24867/IJIEM-2020-1-251

[5] Salman, S. & Alaswad, S. (2021). Designing Reduced Congestion Road Networks via an Elitist Adaptive Chemical Reaction Optimization. Computers & Industrial Engineering, Accepted on Oct 31, 2021. https://doi.org/10.1016/j.cie.2021.107788
