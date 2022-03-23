#!/usr/bin/env python
# coding: utf-8

# ICE evoltion


import pandas as pd
import random
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import time
deltat= 0.1
rD = 0.15
rE = 0.4
rS = 0.01
rF = 0.005
rD = rD*deltat
rE = rE*deltat
rS = rS*deltat
rF = rF*deltat
Ntime = 2000
timeiteration = int(Ntime/deltat)
Re_enter = pd.DataFrame()
Congplasmid = pd.DataFrame()
Rep = pd.DataFrame()
allsum = pd.DataFrame()
Degraded = pd.DataFrame()
cong = pd.DataFrame()
reco = pd.DataFrame()
reg = pd.DataFrame()
asc = pd.DataFrame()
GENEMEAN = pd.DataFrame()
GenVSice = pd.DataFrame()
GenVexciable = pd.DataFrame() # to keep track of ice capable of conjugation
GenVnonxeciable = pd.DataFrame() # to kepp tyrack of ICE uncapable fo conjugation
GenVcongplasmid = pd.DataFrame() # to keep track of conjugative plasmids
rows, cols = (24, 100)
ICE = [[1 for i in range(rows)] for j in range(cols)]
df = pd.DataFrame(ICE)
orf_map = {'ORF-{0}'.format(i+1): 1 for i in range(24)}
df1 = pd.DataFrame([orf_map])
head = df1.columns
df.columns=head
ICEcomplete = pd.DataFrame(np.array([1]*24))
ICEcomplete = ICEcomplete.transpose()
ICEcomplete.columns = head
reco_cols = ['ORF-1','ORF-2','ORF-3','ORF-4']
regu_cols = ['ORF-5','ORF-6','ORF-7','ORF-8','ORF-9','ORF-10','ORF-11','ORF-12']
assc_cols = ['ORF-11']
cong_cols = ['ORF-13','ORF-14','ORF-15','ORF-16','ORF-17','ORF-18','ORF-19','ORF-20','ORF-21','ORF-22','ORF-23','ORF-24']
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
degraded = 0
for t in tqdm(range(timeiteration)):  # time loop
    df = df.reset_index(drop=True)
    conjable = 0
    nonconjable = 0
    for i in range(len(df)-1):
        for j in range(rows):
            rand = random.uniform(0, 1)
            if rD>rand:
                df.iloc[i,j] = 0 
  
    allsum = df.sum(axis = 1)
    congsum = df[cong_cols].sum(axis = 1)
    recosum = df[reco_cols].sum(axis = 1)
    regsum = df[regu_cols].sum(axis = 1)
    conjable = 0
    #nonconjable = 0
    for i in range(len(df)-1):
        if congsum.iloc[i] > 5:
            conjable += 1
    exable = pd.DataFrame([{'Gen':t, 'Number':conjable}])
    GenVexciable = GenVexciable.append(exable)  # excis able
    ###############################################################
    nonconjable = len(df)-conjable
    nonexciable = pd.DataFrame([{'Gen':t, 'Number':nonconjable}])
    GenVnonxeciable = GenVnonxeciable.append(nonexciable) 
    #################################################################
    re_enter = [] # keep track of ICE that can renter into genome
    congplasmid =[] # keep track oof ICE not capable of entering into genome
    for i in range (len(congsum)-1):
        rand = random.uniform(0, 1)
        if congsum.iloc[i] > 5 and recosum.iloc[i] >1 and rand < rE:
            re=df.iloc[i].copy()
            re_enter.append(re)
            #print(congsum.iloc[i], rand, rE, i)
        elif congsum.iloc[i] > 5 and recosum.iloc[i] ==0 and rand < rE:
            plas = df.iloc[i].copy()
            congplasmid.append(plas)
            #print(recosum.iloc[i], rand, rE, i)
            df.iloc[i] = 0
    re_enter = pd.DataFrame(re_enter)
    df = df.append(re_enter) #append to the total population ICEs that were capable to reenter into the genome
    congplasmid = pd.DataFrame(congplasmid)
    allsum = df.sum(axis = 1)
    for i in range (len(allsum)):
        if allsum.iloc[i] ==0:
            df.drop([i], inplace=True)
            degraded += 1 # keep track of in sll generations
    degr = pd.DataFrame([{'Gen':t, 'Number':degraded}])
    Degraded = Degraded.append(degr)  # Number of degraded ICEs 
    ############################################
    rep = []
    asc = df.iloc[:,10] # tetM gene
    for i in range (len(asc)):
        rand = random.uniform(0, 1)
        if rand < rS*asc.iloc[i]:
            newrep = df.iloc[i].copy()
            rep.append(newrep)
    rep = pd.DataFrame(rep)
    df = df.append(rep)
    Re_enter = Re_enter.append(re_enter)
    Congplasmid = Congplasmid.append(congplasmid)
    Rep = Rep.append(rep)
    for i in range(len(df)-1):
        rand = random.uniform(0,1)
        if rand < rF:
            df = df.append(ICEcomplete) 
    step = t*deltat
    if step%5 == 0:
        gmean = df.mean(axis=0)
        gmean = pd.DataFrame(gmean)
        gmean = gmean.transpose()
        gmean.insert(loc=0, column='Time', value=step)
        GENEMEAN = GENEMEAN.append(gmean)
        
    genvsice = pd.DataFrame([{"Gen": t, "Copies": len(df)}])
    GenVSice = GenVSice.append(genvsice)
    #break
    time.sleep(0.1)
GENEMEAN.to_csv(r'genemean.csv',index = False)
df.to_csv(r'ICEPresAbsTry.csv', index = False)
GenVSice.to_csv(r'GenVSice.csv', index = False)
GenVexciable.to_csv(r'GenVexciable.csv',index = False) # ICEs capale of conjugation
Degraded.to_csv(r'resultt/Degraded.csv',index = False) # Degraded ICE
summean = pd.DataFrame()
summean['reco'] = GENEMEAN[reco_cols].sum(axis=1)/4
summean['regu'] = GENEMEAN[regu_cols].sum(axis=1)/8
summean['assc'] = GENEMEAN[assc_cols]/1
summean['cong'] = GENEMEAN[cong_cols].sum(axis=1)/12
plt.figure(figsize=(10, 5),dpi =100)
plt.subplot(2,1,1)
plt.plot(GENEMEAN['Time'], summean['reco'], color = 'blue')
plt.plot(GENEMEAN['Time'], summean['regu'], color = 'green')
plt.plot(GENEMEAN['Time'], summean['assc'], color = 'red')
plt.plot(GENEMEAN['Time'], summean['cong'], color = 'yellow')
plt.legend(['Recombination', 'Regulation', 'Asscesory',  'Conugation'])
plt.xlabel('Time')
plt.ylabel('Avergae number of genes')
plt.title('Pinit ='+str(cols)+',rD ='+str(rD)+',rE = '+str(rE)+',rS='+str(rS)+',rF='+str(rF))
plt.subplot(2, 1, 2)
plt.plot(GenVSice['Gen'], GenVSice['Copies'],  color = 'blue')
plt.plot(GenVexciable['Gen'], GenVexciable['Number'],  color = 'red')
plt.plot(GenVnonxeciable['Gen'], GenVnonxeciable['Number'],  color = 'green')
plt.plot(Degraded['Gen'], Degraded['Number'],  color = 'yellow')
plt.legend(['# of ICE', '# of excisable', '# of non_excisable',  '# of fully degraded ICE'])
plt.xlabel('Generation')
plt.ylabel('Number')
plt.savefig('figure.pdf')
plt.show()
