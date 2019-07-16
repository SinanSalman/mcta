# Markov Chain Traffic Assignment (MCTA)
This is an implementation in python of the MCTA approach proposed by [1] and expanded upon in [2,3,4].

The following table describes the use of each mcta module. Examples of its use can be found in the jupyter files `test_mcta.ipynb` and `test_mcta_runtime.ipynb`



File|Description
---|---
mcta_cli.py|MCTA Command Line Interface (CLI). *This is where you should start.*
mcta.py|MCTA module. *This is where the main MCTA logic reside.*
mcta_edit.py|MCTA network modification module
mcta_rw.py|MCTA input/output data read/write module
mcta_vis.py|MCTA visualization module

# Gnetic Algorithm
`ga.py` and `ga_batch_run.py` implements a genetic algorithm metaheuristic on top of MCTA. resutls are graphed using `ga_polts.ipynb`.

# emissions.py
A module for estimating emissions based on road network traffic conditions. Estimation methedology is based on: Bureau of Public Roads (BPR) curve, COPERT_v5, and TRANSYT7f.

# Referances
[1] Crisostomi, E., Kirkland, S., & Shorten, R. (2011). A Google-like model of road network dynamics and its application to regulation and control. International Journal of Control, 84(3), 633–651. https://doi.org/10.1080/00207179.2011.568005

[2] Salman, S., & Alaswad, S. (2017). Urban road network crisis response management: Time-sensitive decision optimization. Proceedings of the 2017 Industrial and Systems Engineering Conference, 1307–1313. Retrieved from Scopus.

[3] Salman, S., & Alaswad, S. (2018). Alleviating road network congestion: Traffic pattern optimization using Markov chain traffic assignment. Computers & Operations Research, 99, 191–205. https://doi.org/10.1016/j.cor.2018.06.015

[4] Salman, S., & Alaswad, S. (2019). Mitigating the Impact of Congestion Minimization on Vehicles’ Emissions in a Transportation Network. Proceedings for th 25th International Joint Conference on Industrial Engineering and Operations Management.