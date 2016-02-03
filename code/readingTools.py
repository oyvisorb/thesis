import re
import sys
import cyvcf
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.SeqUtils import nt_search
from Bio.Alphabet import generic_dna

mutationTypes = {}

variantVEP =  {}

tripletData = {
    "G[G>A]G/C[C>T]C":0,
    "G[G>C]G/C[C>G]C":0,
    "G[G>T]G/C[C>A]C":0,
    "C[C>G]C/G[G>C]G":0,
    "C[C>A]C/G[G>T]G":0,
    "C[C>T]C/G[G>A]G":0,
    "T[T>G]T/A[A>C]A":0,
    "T[T>C]T/A[A>G]A":0,
    "T[T>A]T/A[A>T]A":0, 
    "A[A>C]A/T[T>G]T":0,
    "A[A>T]A/T[T>A]T":0,
    "A[A>G]A/T[T>C]T":0
}

tripletData2 = {
    "A[A>C]A/T[T>G]T":0,
    "A[A>G]A/T[T>C]T":0,
    "A[A>T]A/T[T>A]T":0,
    "C[C>A]C/G[G>T]G":0,
    "C[C>T]C/G[G>A]G":0,
    "C[C>G]C/G[G>C]G":0,
}


def isRepetitive(pos, ref):
    return True if ref.seq[pos].isupper() else False
     
def storeMutType(mutType):

    if mutType not in mutationTypes.keys():
        mutationTypes.update({mutType:1})
        print "new type added", mutType
    else:
        mutationTypes[mutType] += 1

def storeVariantVEP(mut):

    if mut not in variantVEP.keys():
        variantVEP.update({mut:1})
        print "new variant added", mut
    else:
        variantVEP[mut] += 1
  
def getMutForm(record):
    
    REF = record.REF
    ALT = record.ALT[0]

    if (REF == 'A' and ALT == 'C'): return 'A>C:T>G'
    if (REF == 'A' and ALT == 'G'): return 'A>G:T>C'
    if (REF == 'A' and ALT == 'T'): return 'A>T:T>A'
    if (REF == 'C' and ALT == 'A'): return 'C>A:G>T'
    if (REF == 'C' and ALT == 'G'): return 'C>G:G>C'
    if (REF == 'C' and ALT == 'T'): return 'C>T:G>A'
    if (REF == 'G' and ALT == 'A'): return 'G>A:C>T'
    if (REF == 'G' and ALT == 'C'): return 'G>C:C>G'
    if (REF == 'G' and ALT == 'T'): return 'G>T:C>A'
    if (REF == 'T' and ALT == 'A'): return 'T>A:A>T'
    if (REF == 'T' and ALT == 'C'): return 'T>C:A>G'
    if (REF == 'T' and ALT == 'G'): return 'T>G:A>C'

    print "Error!"
    sys.exit
    
def storeTripletMutType(ref, pos, record, mut):

    if ref.seq[pos-1:pos+2].upper() != ref.seq[pos].upper()*3:

        mut =  mut.split(':')

        if record.REF == 'A':
            tripletData['A['+mut[0]+']A/T['+mut[1]+']T'] += 1
            return True
        if record.REF == 'C':
            tripletData['C['+mut[0]+']C/G['+mut[1]+']G'] += 1
            return True
        if record.REF == 'G':
            tripletData['G['+mut[0]+']G/C['+mut[1]+']C'] += 1
            return True
        if record.REF == 'T':
            tripletData['T['+mut[0]+']T/A['+mut[1]+']A'] += 1
            return True

        print "Error! at store tripletMutType"
        sys.exit()

    return False


def storeTripletMutType2(ref, pos, record, mut):

    if ref.seq[pos-1:pos+2].upper() != ref.seq[pos].upper()*3:

        mut =  mut.split(':')

        if mut[0] == 'A>C' or mut[0] == 'T>G':
            tripletData2['A[A>C]A/T[T>G]T'] += 1
            return True
        if mut[0] == 'A>G' or mut[0] == 'T>C':
            tripletData2['A[A>G]A/T[T>C]T'] += 1
            return True
        if mut[0] == 'A>T' or mut[0] == 'T>A':
            tripletData2['A[A>T]A/T[T>A]T'] += 1
            return True
        if mut[0] == 'C>A' or mut[0] == 'G>T':
            tripletData2['C[C>A]C/G[G>T]G'] += 1
            return True
        if mut[0] == 'C>G' or mut[0] == 'G>C':
            tripletData2['C[C>G]C/G[G>C]G'] += 1
            return True
        if (mut[0] == "C>T") or (mut[0] == 'G>A'):
            tripletData2['C[C>T]C/G[G>A]G'] += 1
            return True

        print "Error! at storeTripletMutType2"
        print mut
        print mut[0] == 'C>T'
        sys.exit()

    return False

def locateTriplet(i, ref):
    
    # check if triplet
    if ref.seq[pos-1:pos+2].upper() != ref.seq[pos].upper()*3:
        return None
    
    # get the whole single base repeat
    pattern = ref.seq[pos].upper()
    j = i
    while re.match(pattern, ref[i].upper()):
        i += 1
    while re.match(pattern, ref[j].upper()):
        j -= 1
    return [j,i]

def getVEPInfo(item, record):
   
    vep = record.INFO['CSQ']
    if vep:
        vepData = vep.split('|')
        for strInfo in filter(None, vepData):
            if re.match(item, strInfo):
                return strInfo
    return None
           
def getVepData(record):

    vep = record.INFO['CSQ']
    vepData = filter(None ,vep.split('|'))
    return vepData
 
def getRecordInfo(item, record):
    
    list1 = []
    for strInfo in record.INFO.keys():
        if re.match(item, strInfo):
            list1.append(record.INFO[strInfo])
    if list1: 
        return list1 
    return None

def locatePatternArea(ref, i, pattern):  
  
    j = i
    while re.match(pattern, ref[i]):
        i += 1
    while re.match(pattern, ref[j]):
        j -= 1

    if i-j > 1:
        return  [j,i]
    else: 
        return None

def locateWord(ref, area, word):

    if area == None: return None
    seq =  ref.seq[area[0]:area[1]]
    words = nt_search(str(seq.upper()), word.upper())
    list1 = []
    for i in words[1:]:
        list1.append([(area[0]+i), (area[0]+i+len(word))])
    if list1: 
        return list1 
    return None

def locateRepWord(ref, i, w):

    j = i
    ch = w[0]
    while re.match(ch, str(ref.seq[i].upper())):
        i+=1
        if ch == w[0]: ch = w[1]
        elif ch == w[1]: ch = w[0]        
    ch = 'C'
    while re.match(ch, str(ref.seq[j].upper())):
        j-=1
        if ch == w[0]: ch = w[1]
        elif ch == w[1]: ch = w[0]
    return [j,i] 


def printContextArea(ref, record, pos, tr, viewl):
    
    print "__"*viewl
    print ref.seq[pos-viewl:pos+viewl]
    if tr:
        print "%s %s" % ((" "*(viewl-(pos-tr[0]))), ref.seq[tr[0]+1:tr[1]])
    print "%s|%s|"% ((" "*(viewl-1)), " ")
    print "%s %s "% ((" "*(viewl-1)), ref.seq[pos])
    print "%s %s "% ((" "*(viewl-1)), record.ALT[0])

####################################################################################

#VCFfile = "testFiles-Genomes1000/chr1.vcf.gz"
#VCFfile = '../vcfFiles/germline/vepped/germline_vep_chr1.vcf.gz'
#VCFfile = "/home/oyzor/Studie/MASTER/testFiles/intersectedFiles/germline_cosmic"
VCFfile = '/home/oyzor/Studie/MASTER/testFiles/intersectedFiles2/intersected_germlineVepped_cosmic_chr1.vcf.gz'

vcf_reader = cyvcf.Reader(open(VCFfile,'r'))
ch_ref = SeqIO.read('../chromFA/chr1.fa', "fasta")


#filename1 = "dict_chromosomes/n_germline_variantVEP_chr1.dict"
#filename2 = "dict_chromosomes/n_germline_mutTypes_chr1.dict"
#filename3 = "dict_chromosomes/n_germline_tripletTypes_chr1.dict"

filename1 = "dict_chromosomes/n_germline_cosmic_variantVEP_chr1.dict"
filename2 = "dict_chromosomes/n_germline_cosmic_mutTypes_chr1.dict"
filename3 = "dict_chromosomes/n_germline_cosmic_tripletTypes_chr1.dict"


nn = 0
for record in vcf_reader:
    
    pos = record.POS-1 
    if (ch_ref.seq[pos].upper() != record.REF):
        print "WTF?"
        sys.exit()

    vepData = getVepData(record)
    storeVariantVEP(vepData[1])
    mut = getMutForm(record)
    storeTripletMutType2(ch_ref, pos, record, mut)
    storeMutType(mut)
    nn+=1

    #printContextArea(ch_ref, record, pos, tr, 25)
    #print "\nChr:%s pos:%s ref:%s alt:%s" % (record.CHROM, record.POS, record.REF, record.ALT[0])
    #print "MUTATION TYPE:", vepData[1]
    #if len(vepData) > 7:
    #    print "CODING TYPE:", vepData[7] 
    #print "REPETITIVE:", isRepetitive(pos, ch_ref)
    #print "TRIPLET++:",ch_ref.seq[tr[0]+1:tr[1]]
    #print "MUTATION FORM:", mut
    #print record.INFO['CSQ']
    #print len(record.INFO['CSQ'].split('|'))


f =  open(filename1,'w')        
for variant in variantVEP.items():
    f.write(str(variant[0])+" "+str(variant[1]))
    f.write('\n')            
f.close()

f2 =  open(filename2,'w')        
for mut in mutationTypes.items():
    f2.write(str(mut[0])+" "+str(mut[1]))
    f2.write('\n')            
f2.close()

f3 =  open(filename3,'w')        
for tri in tripletData2.items():
    f3.write(str(tri[0])+" "+str(tri[1]))
    f3.write('\n')            
f3.close()
print "amount of samples:",nn
