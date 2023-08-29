from demo import *

from pywebio.input import input, FLOAT, input_group
from pywebio.output import put_text, put_processbar, set_processbar, put_file
from pywebio import start_server
import os
from datetime import datetime
import pandas as pd

# #main
# #Deinococcus demo
# ntime=str(datetime.now()).replace(' ','_')
# taxonomy='Deinococcus'
# select_specie='000020685' #
# proteins='WP_012695046.1,WP_012694483.1,WP_049760561.1,WP_012692255.1,WP_012693173.1,WP_162485426.1,WP_012694690.1,WP_012693159.1,WP_012692693.1,WP_012691930.1,WP_012693877.1,WP_012693497.1'
# #DdrO,DdrC,UvrA,UvrB,UvrC,PolA,RecQ,RecQ,RecJ,SSB,UvsE,DciA
# protein_list=proteins.split(',')


#Mora demo
#taxonomy='Acinetobacter'
#or
#taxonomy='Acinetobacter'
#proteins='WP_004925639.1,WP_004928485.1,WP_004920287.1,WP_004923497.1,WP_004922541.1,WP_004920289.1,WP_004930057.1,WP_004923206.1,WP_004930164.1,WP_004927788.1,WP_004923285.1,WP_004930164.1,WP_004923483.1,WP_004929099.1,WP_011182566.1'
#select_specie='000046845'

def web_work():
    put_text('''
    Welcome to Novel cis-element discovery pipeline.
    You can get familiar with our pipeline and reproduce our experiment with setting below：
    
    #Mora demo
    Input Taxonomy:
    Moraxellales or Acinetobacter
    
    NCBI Protein ID:
    WP_004925639.1,WP_004928485.1,WP_004920287.1,WP_004923497.1,WP_004922541.1,
    WP_004920289.1,WP_004930057.1,WP_004923206.1,WP_004930164.1,WP_004927788.1,
    WP_004923285.1,WP_004930164.1,WP_004923483.1,WP_004929099.1,WP_011182566.1
    
    Input Representative Species Refseq ID: 
    GCF_000046845.1
    
    Hope you like it!
    ''')
    info = input_group("Novel cis-element discovery",[
  input('Input Taxonomy', name='tax'),
  input('Input NCBI Protein ID', name='protein_id'),
  input('Input Representative Species', name='select_specie'
    )])
    
    put_text(info['tax'],info['protein_id'],info['select_specie'])
    
    #basic information
    ntime=str(datetime.now()).replace(' ','_')
    taxonomy=info['tax'].replace(' ','')
    select_specie=info['select_specie'].replace(' ','').split('_')[1].split('.')[0]#
    proteins=info['protein_id'].replace(' ','')
    protein_list=proteins.split(',')
    
    if check_before_task(taxonomy,protein_list)==False:
        put_text('The job is rejected. Please check number of taxonomy or proteins')
        raise
    put_text('Check passed!')
    put_processbar('bar')
    set_processbar('bar', 0.5 / 10)
    
    #processing
    if taxonomy not in os.listdir('/data/songshuang/server_find_novel_regulator/0_gbk'):
        download_gbk(taxonomy)
        extract_protein(taxonomy)
        diamond_task(taxonomy)
    set_processbar('bar', 2 / 10)
    get_selected_protein_diamond(select_specie,protein_list,taxonomy,ntime)
    set_processbar('bar', 4 / 10)
    combine_meme_work(ntime)
    set_processbar('bar', 9 / 10)
    motif_for_R(ntime)
    set_processbar('bar', 10 / 10)
    
    #download the result
    os.chdir('/data/songshuang/server_find_novel_regulator/99_User_project/%s/3_meme_result'%ntime)
    os.system('tar -czvf %s.tar.gz *'%(taxonomy))
    content = open('./%s.tar.gz'%taxonomy, 'rb').read() 
    put_file('%s.tar.gz'%taxonomy, content, 'download me')


if __name__ == '__main__':
    start_server(web_work,port=8848)













