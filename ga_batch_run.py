#!/usr/bin/env python

import ga
import json

repetitions = 10
maxtime = 600
results = []

exp = [{'id': 1, 'settings':{'Psz':125, 'Cpr':0.4, 'Mpr':0.12, 'Ktp':0.3}},
       {'id': 2, 'settings':{'Psz':275, 'Cpr':0.4, 'Mpr':0.12, 'Ktp':0.3}},
       {'id': 3, 'settings':{'Psz':125, 'Cpr':0.6, 'Mpr':0.12, 'Ktp':0.3}},
       {'id': 4, 'settings':{'Psz':275, 'Cpr':0.6, 'Mpr':0.12, 'Ktp':0.3}},
       {'id': 5, 'settings':{'Psz':125, 'Cpr':0.4, 'Mpr':0.24, 'Ktp':0.3}},
       {'id': 6, 'settings':{'Psz':275, 'Cpr':0.4, 'Mpr':0.24, 'Ktp':0.3}},
       {'id': 7, 'settings':{'Psz':125, 'Cpr':0.6, 'Mpr':0.24, 'Ktp':0.3}},
       {'id': 8, 'settings':{'Psz':275, 'Cpr':0.6, 'Mpr':0.24, 'Ktp':0.3}},
       {'id': 9, 'settings':{'Psz':125, 'Cpr':0.4, 'Mpr':0.12, 'Ktp':0.7}},
       {'id':10, 'settings':{'Psz':275, 'Cpr':0.4, 'Mpr':0.12, 'Ktp':0.7}},
       {'id':11, 'settings':{'Psz':125, 'Cpr':0.6, 'Mpr':0.12, 'Ktp':0.7}},
       {'id':12, 'settings':{'Psz':275, 'Cpr':0.6, 'Mpr':0.12, 'Ktp':0.7}},
       {'id':13, 'settings':{'Psz':125, 'Cpr':0.4, 'Mpr':0.24, 'Ktp':0.7}},
       {'id':14, 'settings':{'Psz':275, 'Cpr':0.4, 'Mpr':0.24, 'Ktp':0.7}},
       {'id':15, 'settings':{'Psz':125, 'Cpr':0.6, 'Mpr':0.24, 'Ktp':0.7}},
       {'id':16, 'settings':{'Psz':275, 'Cpr':0.6, 'Mpr':0.24, 'Ktp':0.7}},
       {'id':17, 'settings':{'Psz': 50, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.5}},
       {'id':18, 'settings':{'Psz':350, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.5}},
       {'id':19, 'settings':{'Psz':200, 'Cpr':0.3, 'Mpr':0.18, 'Ktp':0.5}},
       {'id':20, 'settings':{'Psz':200, 'Cpr':0.7, 'Mpr':0.18, 'Ktp':0.5}},
       {'id':21, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.06, 'Ktp':0.5}},
       {'id':22, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.30, 'Ktp':0.5}},
       {'id':23, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.1}},
       {'id':24, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.9}},
       {'id':25, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.5}},
       {'id':26, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.5}},
       {'id':27, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.5}},
       {'id':28, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.5}},
       {'id':29, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.5}},
       {'id':30, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.5}},
       {'id':31, 'settings':{'Psz':200, 'Cpr':0.5, 'Mpr':0.18, 'Ktp':0.5}}]

for x in exp:
    for i in range(repetitions):
        print(str(x),'repetition:',i+1)
        r = ga.ga_batchrun(MaxTime=maxtime, nGn=999999, Psz=x['settings']['Psz'], Cpr=x['settings']['Cpr'], Mpr=x['settings']['Mpr'], Ktp=x['settings']['Ktp'], gMpr=0.1)
        print('      fitness=',r['fitness'],' obj_evals=',r['obj_evals'],' runtime=',r['runtime'])
        results.append({'id':x['id'], 'settings':x['settings'], 'repetition':i+1, 'runtime':r['runtime'], 'fitness':r['fitness'], 'obj_evals':r['obj_evals']})

# write results to text file
print("Saving results to 'ga_batch_run.json' file...")
with open('ga_batch_run.json','w') as outputfile:
    json.dump(results,outputfile,indent=3,sort_keys=True)
