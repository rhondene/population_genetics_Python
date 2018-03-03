#Combined HW1 and HW2
#Author: Rhondene J Wint Feb 2018
"""  Read a chr21.vcf from the 1000genomes dataset to record the number of bi-allelic variants per population, number of
segregating sites for a population and the number of global singletons per population. 
This script is written to run on cluster command-line since chr21.vcf file is 11GB  
#on the UNIX terminal  type command: zcat chr21.vcf.gz file | python variant.py geogroup.txt """

import pandas as pd
import os
import sys


#read in population group for each sample from panel file that labels the the Geographic groups for all 2500 indviduals
geotable= pd.read_table(sys.argv[1], header=0, usecols = ['sample','pop'])
geogroup = geotable['pop'].unique() #stores all the populations


inds = geotable['sample'].unique()
ind_dict = dict()
# individual:count dictionary
for ind in inds:
    ind_dict[ind]=0
    
ind_single[ind]= {{k:0 for (k,v) ind_dict.items()}  #counts all individuals who are global singletons

#generate a dictionary of population:individuals key:value pair
pop_dict = dict()
for group in geogroup:   
    sample = geotable[geotable['pop']==group]['sample'] 
    pop_dict[group] = sample  


pop_dict['variants']= {k:0 for (k,v) in pop_dict.items()}  #records total variants for each population; could be useful later 
pop_dict['seg_sites'] = {k:0 for (k,v) in pop_dict.items()}  #count number of segregating sites for eac population
pop_dict['site_tracker'] = {k:0 for (k,v) in pop_dict.items()} #tracks if the current chromosome position was already recorded as a singleton for that population

""" Segregating sites are positions which show differences (polymorphisms) between related genes in a sequence alignment 
That is any sites where at least one individual is a variant.
     E.g. if of the 500 variants are found to belong to ACB population, then ACB is only updated by +1 for that chromosome position.
if pop_dict['site_tracker'][pop]==0 then no seg site has been updated for that population as yet
if pop_dict['site_tracker'][pop]==1 then that chr position has already been recorded as a segregating site"""


#search for population of an individual
def reverse_lookup(dic, value):
    for key in dic:
        if value in dic[key].values:
            return key

#read vcf line by line and record segregating sites per population
n=os.open("gunzip -c chr21.vcf.gz | grep -c '##').readline()

 
    i =0;
    for line in sys.stdin:
        i+=1
        if (i == n):  #this is the line that stores all column headings with IDs of sequenced individuals
            col_ID = sys.stdin.readline().split('\t') 
            col_ID = col_ID[9:]   #only want columns  with IDs of individuals
        if i > n:
            row = sys.stdin.readline().split('\t')
            try:  #i got an index error for some reason when I ran the whole vcf file, but did not appear on subset of data
                
                if len(row[3])==1 and len(row[4])==1: #bi-allelic sites, a single allele in ref column and alt column 
                    row=row[9:]
                #reset trackers to zero before new row  
                    pop_dict['site_tracker'] = {k:0 for (k,v) in pop_dict['site_tracker'].items()} 
                    single_list=[]
                    for x in range(len(row)):  
                        if (row[x]=='0|1' or row[x]=='1|0' or row[x]=='1|1'):
                            ind_dict[col_ID[x]] +=1
                        
                            pop = reverse_lookup(pop_dict,col_ID[x])  #returns population key of indivdiual
                            pop_dict['variants'][pop]+=1
                        #update segregating site if not yet done
                            if (pop_dict['site_tracker'][pop]==0 ):
                                pop_dict['seg_sites'][pop]+=1
                                pop_dict['site_tracker'][pop]+=1 #update tracker for that population
                  
                        if (row[x]=='0|1' or row[x]=='1|0'):
                            if (len(single_list)<2):
                                single_list.append(col_ID[x])  #stores singleton's ID
                    #update the individual corresponding to the singleton
                    if (len(single_list)==1):
                        ind_single[single_list[0]]+=1 
                                  
            except IndexError:   #restart loop
                continue                         
       # if (i == 2000): break; for testing if code works on subset of data     
   
#could've converted the dicts directily to csv but want to avoid extra steps of fixing indexing issues
seg_sites = pd.DataFrame(list(pop_dict['seg_sites'].items()), columns=['Population', 'Seg_Sites'])
var_sites = pd.DataFrame(list(ind_dict.items()), columns=['Individual', 'Var_Sites'])
variants = pd.DataFrame(list(pop_dict['variants'].items()), columns = ['Population', 'No_Variant'])
singletons = pd.DataFrame(list(ind_single.items()), columns = ['Population', 'No_Singleton'])
                  
#map individual variants to their population
var_sites.sort_values(var_sites.columns[0], inplace=True)
geotable.sort_values(geotable.columns[0], inplace=True)
var_sites['Population'] = geotable['pop']
                 
singletons.sort_values(singletons.columns[0], inplace=True)
singletons['Population'] = geotable['pop']
#Cluster outputs csv file to access later for visualisation in Jupyter notebook

seg_sites.to_csv('seg_sites.csv', sep=',', index=False)
singletons.to_csv('singletons.csv', sep=',', index =False)
var_sites.to_csv('var_sites.csv', sep=',', index=False)
variants.to_csv('variants.csv', sep=',', index=False)

      """" Plotting Figure 1b from 1000Genomes Paper in Jupyter Notebook"""
import matplotlib.pyplot as plt
import matplotlib

%matplotlib notebook
matplotlib.get_backend()
          
import seaborn as sns
         """ Plot Number of Variants on Chromosome 21 Per Population"""
var_sites2= pd.read_csv('var_sites2.csv', header=0, delimiter =',')
var_sites2.tail(5)
          
plt.figure(figsize=(10,6))
sns.stripplot(x='Population', y='Var_Sites', data= var_sites2, size=5)  #store an array of all the strip objects

# Modify Graph Aesthetics
sns.set_style('white')
sns.set_palette("Paired", 4)
sns.despine()
plt.style.use('seaborn-colorblind')

ax = plt.gca()
ax.set_xlabel("Individual")
ax.set_ylabel("Variant Sites Per Genome")
ax.set_title(" Number of Variant Sites on Chr.21 In Individuals from Each Populations ")
plt.xticks(rotation = 'vertical')

#Plot ingoring the individual  who had 0 var_sites
plt.figure(figsize=(10,6))
sns.stripplot(x='Population', y='Var_Sites', data= var_sites2[var_sites2.Var_Sites !=0], size=5)  #store an array of all the strip objects
sns.despine()
# Modify Graph Aesthetics
sns.set_style('white')
sns.set_palette("Paired", 4)
plt.style.use('seaborn-colorblind')

ax = plt.gca()
ax.set_xlabel("Individual")
ax.set_ylabel("Variant Sites Per Genome")
ax.set_title(" Number of Variant Sites on Chr.21 In Individuals from Each Populations ")
plt.xticks(rotation = 'vertical')
          
          
 """Plot Number of Global Singletons for Each Population
Plot Number of Segregating Sites for Each Population"""
          
 seg_sites = pd.read_csv('seg_sites2.csv', header=0, delimiter=',')
singletons = pd.read_csv('singletons.csv', header=0, delimiter=',')
 variants = pd.read_csv('variants2.csv', header=0, delimiter=',')
  plt.figure(figsize=(10,6))
strips =sns.stripplot(x='Population', y='No_Singleton', data= singletons, size=5)  #store an array of all the strip objects

# Modify Graph Aesthetics
sns.set_style('white')
sns.set_palette("Paired", 4)
plt.style.use('seaborn-colorblind')
plt.xticks(rotation = 'vertical')
plt.ylabel('Number of Global Singletons')
plt.xlabel('Individual')
plt.title('Number of Global Singles in Each Population')
 
          
      #Visualise and Compare How Singletons, Variants and Segregating Sites Vary Among Sequenced Population 
plt.figure(figsize=(10,10))

ax0 = plt.subplot2grid((8,3), (0,0), rowspan=3, colspan=3)                       
ax0=sns.barplot(x='Population', y='Seg_Sites', data= seg_sites)  #store an array of all the strips
ax0.set_ylabel("Segregating Sites Per Genome")
ax0.set_title(" Number of Segregating Sites on Chr.21 for Each Population")
for tick in ax0.get_xticklabels():
        tick.set_rotation('vertical')
      

ax1 = plt.subplot2grid((8,3), (3,0), rowspan=3, colspan=3)
ax1 = sns.barplot(x= 'Population', y='No_Singleton', data= singletons )
ax1.set_ylabel("Singletons Per Genome")
ax1.set_title(" Number of Global Singletons on Chr.21 For Each Populations")
for tick in ax1.get_xticklabels():
        tick.set_rotation('vertical')

ax2 = plt.subplot2grid((8,3), (6,0), rowspan=4, colspan=3) #starts at the second row and takes up 1st column to rest of column space 
ax2 = sns.barplot (x= 'Population', y='No_Variant', data= variants )
ax2.set_ylabel("Total Variants Per Genome")
ax2.set_title(" Total Variants on Chr.21 For Each Populations")
for tick in ax2.get_xticklabels():
        tick.set_rotation('vertical')


sns.set_style('white')
sns.set_palette("Paired", 6)
plt.style.use('seaborn-dark-palette')

plt.tight_layout()

plt.show()

          
          
