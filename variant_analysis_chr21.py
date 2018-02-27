#Combined HW1 and HW2
#Author: Rhondene J Wint
"""  Read a chr21.vcf from the 1000genomes dataset to record the number of bi-allelic variants per population, number of
segregating sites for a population and the number of global singletons per population"""

import pandas as pd
import os
import sys


#read in population group for each sample from panel file
geotable= pd.read_table(sys.argv[1], header=0, usecols = ['sample','pop'])
geogroup = geotable['pop'].unique() #stores all the populations


inds = geotable['sample'].unique()
ind_dict = dict()
# individual:count dictionary
for ind in inds:
    ind_dict[ind]=0

#generate a dictionary of population:individuals key:value pair
pop_dict = dict()
for group in geogroup:   
    sample = geotable[geotable['pop']==group]['sample'] 
    pop_dict[group] = sample  


pop_dict['variants']= {k:0 for (k,v) in pop_dict.items()}  #records total variants for each population; could be useful later 
pop_dict['seg_sites'] = {k:0 for (k,v) in pop_dict.items()}  #count number of segregating sites for eac population
pop_dict['site_tracker'] = {k:0 for (k,v) in pop_dict.items()} #tracks if the current chromosome position was already recorded as a singleton for that population
pop_dict['singleton'] = {k:0 for (k,v) in pop_dict.items()} # count the number of global singletons that belong to a population 

#search for population of a individual
def reverse_lookup(dic, value):
    for key in dic:
        if value in dic[key].values:
            return key

#read vcf line by line and record segregating sites per population
n=252  # header lines to skip; n = gunzip -c chr21.vcf.gz | grep -c '##' or n =os.open("gunzip -c chr21.vcf.gz | grep -c '##').readline()

with open(sys.argv[2], 'r+') as fh:
    i =0;
    for line in fh:
        i+=1
        if (i == n):  
            col_ID = fh.readline().split('\t') 
            col_ID = col_ID[9:]
        if i > n:
            row = fh.readline().split('\t')
            try:  #i got an index error for some reason when I ran the whole vcf file, but did not appear on subset of data
                
                if len(row[3])==1 and len(row[4])==1: #bi-allelic sites, a single allele in ref column and alt column 
                    row=row[9:]
                #reset trackers to zero before new row  
                    pop_dict['site_tracker'] = {k:0 for (k,v) in pop_dict['site_tracker'].items()} 
                    single_list=[]
                    for x in range(len(row)):  #could also do for col in row
                        if (row[x]=='0|1' or row[x]=='1|0' or row[x]=='1|1'):
                            ind_dict[col_ID[x]] +=1
                        
                            pop = reverse_lookup(pop_dict,col_ID[x])  #returns population key of indivdiual
                            pop_dict['variants'][pop]+=1
                        #update segregating site if not yet done
                            if (pop_dict['site_tracker'][pop]==0 ):
                                pop_dict['seg_sites'][pop]+=1
                                pop_dict['site_tracker'][pop]+=1 #update tracker for that population
                        #determines if global singleton 
                            if (len(single_list)<2):
                                single_list.append(col_ID[x])  #stores singleton's ID
                    #and map the singleton to its population
                    if (len(single_list)==1):
                        pop = reverse_lookup(pop_dict, single_list[0]) #looks up population for the indivdiual
                        pop_dict['singleton'][pop]+=1
                                  
            except IndexError:   #restart loop
                continue                         
        if (i == 2000): break;     
   
#could've converted the dicts directily to csv but want to avoid indexing issues
seg_sites = pd.DataFrame(list(pop_dict['seg_sites'].items()), columns=['Population', 'Seg_Sites'])
singletons = pd.DataFrame(list(pop_dict['singleton'].items()), columns = ['Population', 'No_Singleton'])
var_sites = pd.DataFrame(list(ind_dict.items()), columns=['Individual', 'Var_Sites'])
variants = pd.DataFrame(list(pop_dict['variants'].items()), columns = ['Population', 'No_Variant'])

#map individual variants to their population
var_sites.sort_values(var_sites.columns[0], inplace=True)
geotable.sort_values(geotable.columns[0], inplace=True)
var_sites['Population'] = geotable['pop']

#store as csv file for later visualisation

seg_sites.to_csv('seg_sites.csv', sep=',', index=False)
singletons.to_csv('singletons.csv', sep=',', index =False)
var_sites.to_csv('var_sites.csv', sep=',', index=False)
variants.to_csv('variants.csv', sep=',', index=False)


