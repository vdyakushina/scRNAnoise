import pandas as pd
import os
from subprocess import run
import sys

ref_homo="/media/CACHE-MUT/starsolo/ref_Homo"
ref_mouse="no"
working_dir=sys.argv[1]
db=pd.read_csv('%s' % working_dir+'sample_annotated.tsv', sep='\t', header=None)
chem={'GSE121521':'v2', 'GSE157829': 'v2', 'GSE128334':'v2', 'GSE145708':'v3', 'GSE137478':'v2', 'GSE154873':'v3', 'GSE169751':'v2', 'GSE165238':'v3', 'GSE160563':'v3', 'GSE150482':'v3', 'GSE135922':'v3', 'GSE103983':'DropSeq', 'GSE97104':'DropSeq', 'GSE96827':'DropSeq', 'GSE93374':'DropSeq'}
db=db[db[7]=="Homo sapiens"]
db=db[~((db[0]=='GSE160563')&(~db[5].isin(['SRR12957805', 'SRR12957806', 'SRR12957807', 'SRR12957808', 'SRR12957809', 'SRR12957810'])))]
db=db[db[0].isin(chem.keys())]

for index, line in db.iterrows():
  experiment=working_dir+line[0]
  sample=working_dir+line[0]+line[5]
  if chem[line[0]]==v2:
    barcodes='/media/CACHE-MUT/starsolo/737K-august-2016.txt'
  elif chem[line[0]]==v3:
    barcodes='/media/CACHE-MUT/starsolo/3M-february-2018.txt'
  else:
    barcodes=None
  for path in line[8].split(';'):
    run("mkdir %s" % sample, shell = True)
    run("wget -P %s %s" % (sample, path), shell = True)
  files=sorted([db.loc[46][8].split(';')[0].split('/')[-1], db.loc[46][8].split(';')[1].split('/')[-1]])
  if "bam" in files[0]:
    module="bam"
  else:
    module="fastq"
  if db[7]=="Homo sapiens":
    ref_dir=ref_homo
  else:
    ref_dir=ref_mouse
  if module=="fastq":
    run("STAR --genomeDir %s --readFilesIn %s %s --readFilesCommand zcat --soloType CB UMI Simple --soloCBwhitelist %s --outFileNamePrefix %s" %
        (ref_dir, sample+"/"+files[1], sample+"/"+files[0], barcodes, sample), shell = True)
  else:
    run("STAR --genomeDir %s --readFilesIn %s --readFilesType SAM SE --soloType CB UMI Simple --soloCBwhitelist %s --outFileNamePrefix %s" %
        (ref_dir, sample+"/"+files[0], barcodes, sample), shell = True)