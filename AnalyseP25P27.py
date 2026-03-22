#!/usr/bin/python3
#-*- coding : utf-8 -*-

__authors__ = ("Lucien Maurau")
__contact__ = ("lucien.maurau@etu.umontpellier.fr")
__version__ = "0.0.1"
__date__ = "11/02/2026"
__licence__ ="This program is free software: you can redistribute it and/or modifyit under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>."

import os;
import sys;

def check_file(file):
    if not os.path.exists(file): #check if the file exists
        print(f"Error : the path '{file}' does not exist.")
        return False
    if not os.path.isfile(file): #check if the path is a file
        print(f"Error : the path '{file}' is not a file.")
        return False
    size = os.path.getsize(file)
    if size == 0: #check if the file is empty
        print(f"Error : the file '{file}' is empty.")
        return False
    print(f"OK: the file '{file}' exists and is not empty (size: {size} bytes)") #confirm the file is valid
    return True

def VcfID(input_file,ext):
    names=["P15NP90Global"+ext+".vcf","P90NP15Global"+ext+".vcf","P25NP27Global"+ext+".vcf","P27NP25Global"+ext+".vcf","P27NP25NP90Global"+ext+".vcf"]
    for s in names:
        VcfWriter(s)
    with open(input_file,"r") as Str:
        S=Str.readlines()
        for Ligne in S:
            if Ligne[0] != '#':
                line=Ligne.split("\t")
                if float(line[5])>30.00:
                    IDs=IDclassifier(line)
                    Wlist=VeenClassifier(IDs)
                    for s in range(0,4):
                        if Wlist[s]:
                            WriteLine(names[s],Ligne)
                    if Wlist[3] & Wlist[4]:
                        WriteLine(names[4],Ligne)
    return 0

def IDclassifier(line):
    ID=[[0 for s in range(5)] for i in range(7)]
    s=9
    Q=0
    q=0
    while s < len(line):
        if line[s][0] != '.':
            ID[q][Q]+=1
        s+=1
        Q+=1
        if Q == 5:
            Q=0
            q+=1
    return ID

def VeenClassifier(ID):
    P15NP90,P90NP15,P25NP27,P27NP25,NP90=False,False,False,False,False
    for s in range(4):
        if ID[0][s]!=0:
            if ID[-1][s]==0:
                P15NP90=True
        else :
            if ID[-1][s]!=0:
                P90NP15=True
    for s in range(4):
        if ID[1][s]!=0:
            if ID[2][s]==0:
                P25NP27=True
        else :
            if ID[2][s]!=0:
                P27NP25=True
    if P27NP25:
        for s in range(4):
            if ID[-1][s]==0:
                NP90=True
    return [P15NP90,P90NP15,P25NP27,P27NP25,NP90]

def VcfWriter(name):
    with open(name,"w") as text:
        text.write("##fileformat=VCFv4.1\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##medaka_version=1.11.3\n##contig=<ID=DQ657948.1,length=272677>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Medaka genotype\">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Medaka genotype quality score\">\n##bcftools_filterVersion=1.21+htslib-1.21\n##bcftools_filterCommand=filter -i QUAL>=30 -o P15-1.trimed1000.snp_filtered.vcf -O v P15-1.trimed1000.snp.vcf; Date=Wed Feb 18 11:12:19 2026\n##bcftools_mergeVersion=1.21+htslib-1.21\n##bcftools_mergeCommand=merge --force-samples -O v -o P15_1-5_merged.vcf P15-1.trimed1000.snp_filtered.vcf.gz P15-2.trimed1000.snp_filtered.vcf.gz P15-3.trimed1000.snp_filtered.vcf.gz P15-4.trimed1000.snp_filtered.vcf.gz P15-5.trimed1000.snp_filtered.vcf.gz; Date=Wed Feb 18 11:59:02 2026\n##bcftools_mergeCommand=merge -O v -o Base_froid.vcf P15_1-5_renamed.vcf.gz P25_1-5_renamed.vcf.gz P27_1-5_renamed.vcf.gz P30_1-5_renamed.vcf.gz P50_1-5_renamed.vcf.gz P65_1-5_renamed.vcf.gz P90_1-5_renamed.vcf.gz; Date=Wed Feb 18 17:03:59 2026\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	P15-1	P15-2	P15-3	P15-4	P15-5	P25-1	P25-2	P25-3	P25-4	P25-5	P27-1	P27-2	P27-3	P27-4	P27-5	P30-1	P30-2	P30-3	P30-4	P30-5	P50-1	P50-2	P50-3	P50-4	P50-5	P65-1	P65-2	P65-3	P65-4	P65-5	P90-1	P90-2	P90-3	P90-4	P90-5\n")
    return 0

def WriteLine(name,line):
    with open(name,"a+") as file:
        file.write(line)

def main(argv):
    if len(argv) != 2:
        print("Erreur, mauvais nombre d'arguments")
        return None
    input_file_1=argv[0]
    ext=argv[1]
    if check_file(input_file_1): 
        VcfID(input_file_1,ext)
    return 0 
        
if __name__ == "__main__":
    main(sys.argv[1:]) 