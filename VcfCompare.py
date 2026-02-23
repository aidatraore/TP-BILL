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

def VcfID(input_file):
    with open(input_file,"r") as Str:
        S=Str.readlines()
        Size=len(S)
        IdList=[0 for i in range(Size)]
        i=0
        for Ligne in S:
            if Ligne[0] != '#':
                IdList[i-6]=Ligne.split("\t")[1]
            i+=1
    return IdList[:-1],S

def VariantDetection(ID1,ID2,l,n):
    Pos=int(l[1])
    Qual=float(l[5])
    if Qual >29.00:
        for s in range(max(Pos-n//2,0),Pos+n//2):
            if str(s) in ID2:
                if str(s) not in ID1:
                    return True
    return False

def VcfWritter(ID1,ID2,S):
    i=0
    with open("Output.vcf","a+") as Output:
        for Ligne in S:
            if Ligne[0] != '#':
                l=Ligne.split("\t")
                if VariantDetection(ID1, ID2, l, 10):
                    Output.write(Ligne)
            i+=1
    return 0

def VcfWriter():
    with open("Output.vcf","w") as text:
        text.write("##fileformat=VCFv4.1\n##medaka_version=1.11.3\n##contig=<ID=DQ657948.1,length=272677>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Medaka genotype\">\n##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Medaka genotype quality score\">\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n")
    return 0

def main(argv):
    if len(argv) != 2:
        print("Erreur, mauvais nombre d'arguments")
        return None
    input_file_1=argv[0]
    input_file_2=argv[1]
    if check_file(input_file_1) & check_file(input_file_2): 
        VcfWriter()
        ID1,S1 = VcfID(input_file_1)
        ID2,S2 = VcfID(input_file_2)
        VcfWritter(ID1,ID2,S2)
        
if __name__ == "__main__":
    main(sys.argv[1:]) 
        
