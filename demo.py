#!/usr/bin/env python
# coding: utf-8


def check_before_task(taxonomy,protein_list):
    os.chdir('/data/songshuang/server_find_novel_regulator')
    if len(protein_list)>15:
        print('Too many proteins, the job is rejected.')
        return False
    #check taxonomy
    os.system('datasets download genome taxon %s --assembly-source refseq --include gbff --reference  --preview > preview'%(taxonomy))
    data=open('preview').read()
    tax_count=data.split(',')[1].split(':')[1]
    print(tax_count)
    os.system('rm preview')
    if int(tax_count) < 5 or int(tax_count) >200 :
        print('Species count is %d, the job is rejected.'%int(tax_count))
        return False
    print('Check passed')




#download from ncbi
import os
def download_gbk(taxonomy):
    os.chdir('/data/songshuang/server_find_novel_regulator/')
    os.system('mkdir /data/songshuang/server_find_novel_regulator/0_gbk/'+taxonomy)
    os.system('mkdir /data/songshuang/server_find_novel_regulator/0_protein/'+taxonomy)
    os.system('mkdir /data/songshuang/server_find_novel_regulator/0_promoter/'+taxonomy)
    
    
    
    os.chdir('/data/songshuang/server_find_novel_regulator/0_gbk/'+taxonomy)
    os.system('datasets download genome taxon %s --assembly-source refseq --include gbff --reference'%(taxonomy))
    os.system('unzip ncbi_dataset.zip')
    os.system('rm -rf ncbi_dataset.zip README.md ncbi_dataset/data/dataset_catalog.json ncbi_dataset/data/assembly_data_report.jsonl')
    
    for folder in os.listdir('ncbi_dataset/data'):
        genome_id=folder.split('.')[0].split('_')[1]
        os.system('mv ncbi_dataset/data/%s/genomic.gbff ./%s'%(folder,genome_id))
        
    os.system('rm -rf ncbi_dataset')
    


#extract protein & promoter
import os
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "shuang_s@zju.edu.cn"
def all(genome_id):
    try:
        protein_data=''
        promoterdata=''
        #print('running')
        for seq_record in SeqIO.parse(genome_id, "genbank"):
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    if 'protein_id' not in feature.qualifiers:
                        continue
                    #print(feature.type)
                    #print(feature.location)
                    protein_seq = feature.qualifiers['translation'][0]

                        
                    if '-' in str(feature.location):
                        promoter_seq = seq_record.seq[(feature.location.end):(feature.location.end + 200)].reverse_complement()
                    else:
                        promoter_seq = seq_record.seq[(feature.location.start - 200):(feature.location.start)]
                        
                    if len(protein_seq) < 30 or len(promoter_seq) < 60:
                        continue
                        
                    try:
                        protein_data += '>%s\n%s\n' % (genome_id+'__'+feature.qualifiers['protein_id'][0] ,protein_seq)#.split('_')[0]
                        promoterdata += '>%s\n%s\n' % (genome_id+'__'+feature.qualifiers['protein_id'][0] ,promoter_seq)#.split('_')[0]

                    except:
                        protein_data += '>%s\n%s\n' % ('unknown' , protein_seq)
                        promoterdata += '>%s\n%s\n' % ('unknown' , promoter_seq)
                        


        f=open('%s.faa'%(genome_id),'w+')
        f.write(protein_data)
        f.close()
        
        f=open('%s.fa'%(genome_id),'w+')
        f.write(promoterdata)
        f.close()
        print('finished')
    except:
        raise
        print(genome_id)



import os
def extract_protein(tax):
    os.chdir('/data/songshuang/server_find_novel_regulator/0_gbk/'+tax)
    from multiprocessing import Pool
    p=Pool(32)
    p.map(all,os.listdir('/data/songshuang/server_find_novel_regulator/0_gbk/'+tax))
    os.system('mv *.fa /data/songshuang/server_find_novel_regulator/0_promoter/'+tax)
    os.system('mv *.faa /data/songshuang/server_find_novel_regulator/0_protein/'+tax)




#build diamond database
def diamond_task(tax):
    os.chdir('/data/songshuang/server_find_novel_regulator/0_protein/'+tax)
    os.system('cat * > all_protein.faa')
    os.system('diamond makedb --in all_protein.faa -p 32 -d all_protein.dmnd')
    



#find promoter sequence
import os
def find_promoter(genome_id,protein_id,tax,form):
    if form=='protein':
        f=open('/data/songshuang/server_find_novel_regulator/0_protein/%s/%s.faa'%(tax,genome_id))
    else:
        f=open('/data/songshuang/server_find_novel_regulator/0_promoter/%s/%s.fa'%(tax,genome_id))
    data=''
    write=False
    for line in f:
        if (line[0]=='>'):
            if write==True:
#                 if len(data)<=50:
#                     print('here!')
#                print(data)
                return data
            if (protein_id in line):
                write=True
        if write==True:
            data+=line
    return '>%s\nAAAAAA\n'%(protein_id)




import os
def get_selected_protein_diamond(select_specie,protein_list,tax,ntime):
    os.system('mkdir /data/songshuang/server_find_novel_regulator/99_User_project/'+ntime)
    os.system('mkdir /data/songshuang/server_find_novel_regulator/99_User_project/%s/0_selected_protein'%ntime)
    os.system('mkdir /data/songshuang/server_find_novel_regulator/99_User_project/%s/1_diamond_result'%ntime)
    os.system('mkdir /data/songshuang/server_find_novel_regulator/99_User_project/%s/1_diamond_result'%ntime)
    os.system('mkdir /data/songshuang/server_find_novel_regulator/99_User_project/%s/2_promoter'%ntime)
    os.system('mkdir /data/songshuang/server_find_novel_regulator/99_User_project/%s/3_meme_result'%ntime)
    

    
    import pandas as pd
    #diamond
    for protein in protein_list:
        os.chdir('/data/songshuang/server_find_novel_regulator/99_User_project/%s/0_selected_protein'%ntime)
        f=open(protein,'w+')
        f.write(find_promoter(select_specie,protein,tax,'protein'))
        f.close()
        os.system('diamond blastp -d /data/songshuang/server_find_novel_regulator/0_protein/%s/all_protein.dmnd -q %s -o ../1_diamond_result/%s -t /dev/shm -b10.0 -k100'%(tax,protein,protein))

    #analysis diamond result
        diamond_result=pd.read_csv('../1_diamond_result/%s'%(protein),sep='\t|;',engine='python',header=None)
        diamond_result=diamond_result.loc[diamond_result.iloc[:,10]<1e-30,:]
    #extract promoter
        f=open('/data/songshuang/server_find_novel_regulator/99_User_project/%s/2_promoter/%s'%(ntime,protein),'w+')
        diamond_protein_list=list(diamond_result.iloc[:,1])
        for diamond_protein in diamond_protein_list:
            genome=diamond_protein.split('__')[0]
            protein_id=diamond_protein.split('__')[1]
            f.write(find_promoter(genome,protein_id,tax,'promoter'))
        f.close()
        
        ##cd-hit the promoter
        os.chdir('/data/songshuang/server_find_novel_regulator/99_User_project/%s/2_promoter'%ntime)
        os.system('cd-hit -i %s -c 0.6 -o cd-hit_%s -T 6 -n 2'%(protein,protein))
        os.system('mv cd-hit_%s %s'%(protein,protein))
        os.system('rm cd-hit*')
        



def meme_work(ntime,com):
    os.chdir('/data/songshuang/server_find_novel_regulator/99_User_project/%s/2_promoter'%(ntime))
    file_name='_'.join(com)
    #print(file_name)
    data=''
    for name in com:
        data+=open(name).read()
    f=open('/dev/shm/%s/%s'%(ntime,file_name),'w+')
    f.write(data)
    f.close()
    
    os.chdir('/dev/shm/%s'%ntime)
    
    
    os.system('/home/songshuang/miniconda3/envs/meme/bin/meme %s -oc /data/songshuang/server_find_novel_regulator/99_User_project/%s/3_meme_result/%s -dna -minw 12 -maxw 25 -mod zoops -evt 1e-40 -p 4 -nmotifs 3'%(file_name,ntime,file_name))
    
    if len(os.listdir('/data/songshuang/server_find_novel_regulator/99_User_project/%s/3_meme_result/%s'%(ntime,file_name))) < 5:
        os.system('rm -r /data/songshuang/server_find_novel_regulator/99_User_project/%s/3_meme_result/%s'%(ntime,file_name))
        print(file_name+' deleted')
    




def combine_meme_work(ntime):
    work=os.listdir('/data/songshuang/server_find_novel_regulator/99_User_project/%s/2_promoter'%(ntime))
    from itertools import combinations
    works=combinations(work,3)
    os.chdir('/data/songshuang/server_find_novel_regulator/99_User_project/%s/2_promoter'%(ntime))
    os.system('mkdir /data/songshuang/server_find_novel_regulator/99_User_project/%s/3_meme_result'%(ntime))
    os.system('mkdir /dev/shm/%s'%ntime)
    from multiprocessing import Pool
    from functools import partial
    p=Pool(10)
    p.map(partial(meme_work, ntime),works)
    os.system('rm -r /dev/shm/%s'%ntime)


def motif_for_R(ntime):
    os.chdir('/data/songshuang/server_find_novel_regulator/99_User_project/%s/3_meme_result'%(ntime))
    os.system('/home/songshuang/.conda/envs/R/bin/Rscript /home/songshuang/code/DdaA文章/99_server_for_novel/motif.R')

