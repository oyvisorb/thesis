import os
import sys
import numpy as np
import gzip
from pybedtools import BedTool, cleanup
#pybedtools.set_tempdir('/usit/abel/u1/oyvisorb/temp')

import clearStages

path = 'stage1/'
path2 = 'stage2/'



def mergeVCFfiles(chrom, A, B, file1, file2):

    ## Write headers to a file. Needed later for creating a correct VCF of the intersected files. 
    header1 = gzip.open((file1+chrom+'.vcf.gz'), 'r')
    header2 = gzip.open((file2+chrom+'.vcf.gz'), 'r')

    headerFile = (path+'HEADER_'+A+'_'+B+'_.vcf.gz')
    f = gzip.open(headerFile, 'ab')

    for line in header1:
        if line[0:2] != '##':
            break;
        f.write(line)
    header1.close()

    for line in header2:
        if line[0:2] != '##':
            break;
        f.write(line)
    header1.close()

    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO');
    f.close()

    ## Intersects files -LOJ

    a = BedTool((file1+chrom+'.vcf.gz'))
    b = BedTool((file2+chrom+'.vcf.gz'))

    a_b = a.intersect(b, header=True,loj=True)
    ab  = a_b.saveas((path+'LOJ_'+A+'_'+B+'_chr'+chrom+'.vcf.gz'))

    print (ab.fn)
    cleanup(remove_all=True)

    ## Intersects files -V

    a = BedTool((file1+chrom+'.vcf.gz'))
    b = BedTool((file2+chrom+'.vcf.gz'))

    b_a = b.intersect(a, header=True, v=True)
    ba  = b_a.saveas((path+'V_'+A+'_'+B+'_chr'+chrom+'.vcf.gz'))

    print (ba.fn)
    cleanup(remove_all=True)


    ## CAT LOJ file and -v File

    LOJ = path+'LOJ_'+A+'_'+B+'_chr'+chrom+'.vcf.gz'
    V = path+'V_'+A+'_'+B+'_chr'+chrom+'.vcf.gz'
    out = path+'concated_'+A+'_'+B+'_chr'+chrom+'.vcf.gz'

    os.system(("vcf-concat "+LOJ+" "+V+" | gzip -c > "+out))
    #print "ok"



    ## correct to vcf, merge equal type samples
    out2 = 'stage2/unsorted_'+A+'_'+B+'_chr'+chrom+'.vcf.gz'

    import correctToVCF
    goVCF = correctToVCF.CorrectToVCF()
    goVCF.writeHeader(headerFile, out2)
    goVCF.correctFile(out, out2)

    ## sort the VCF file

    os.system(("vcf-sort "+out2+"| gzip -c > "+("germlineVCF/"+A+'_'+B+'_chr'+chrom+'.vcf.gz')))

    cleanup(remove_all=True)

chromosomes = ['1','2','3','4','5','6','7','8','9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

for i in chromosomes[1::]:

    chrom = str(i)

    A = 'dbSNP'
    B = 'Genomes1000'
    C = 'ExAC'

    #file1 = 'testFiles-dbSNP/chr'
    #file2 = 'testFiles-Genomes1000/chr'
    #file3 = 'testFiles-ExAC/chr'

    file1 = '../vcfFiles/dbSNP/chromosomes-dbSNP/chr'
    file2 = '../vcfFiles/1000genomes/chromosomes-1000genomes/chr'
    file3 = '../vcfFiles/ExAC/chromosomes-ExAC/chr'

    mergeVCFfiles(chrom, A, B, file1, file2)

    AB = 'dbSNP_Genomes1000'
    file4 = 'germlineVCF/'+A+'_'+B+'_chr'

    os.system("python clearStages.py")
    
    mergeVCFfiles(chrom ,C, AB, file3, file4)
    
    os.system("python clearStages.py")

print "Succes!"
