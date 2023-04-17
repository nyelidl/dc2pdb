#!/usr/bin/python3

# created by Duan 2022.10.14
import pandas as pd
import numpy as np
import argparse
import os


#parser.add_argument('-l', '--location', type=str, help='ENTER YOUR LOCATION', required=True)

#os.system("awk '{print $2}' " + args.location + " > arg_atomNo.dat")
print("---------------------------------------------------------")
print("Welcome to dcfactor! \ncreated by Duan 2022.10.14")
print("change temperature factor to DC total energy to PDB file")
print("---------------------------------------------------------")
print("Residue Name in PDB need starting from 1!")
print("need package: pandas, numpy. \nyou can install package by	   \n \npip3 install numpy  \npip3 install pandas  \n")
parser = argparse.ArgumentParser(description=" ")

print("if your pdb files don't has chain please add by:		  \n\nsed -i '1,$s/^\(ATOM.\{17\}\)/\\1A/' ***.pdb        \nsed 's/ A / A/g' ***.pdb 	 \n	")
parser.add_argument('-p', '--pdbfile', type=str, help='ENTER YOUR PDB FILE LOCATION', required=True)
parser.add_argument('-d', '--dcfile', type=str, help='ENTER YOUR DC FILE LOCATION', required=True)
parser.add_argument('-o', '--output', type=str, help='ENTER YOUR OUTPUT LOCATION', default="dl.pdb")
args = parser.parse_args()

os.system("sed '/^Sidechain/,$d' " + args.dcfile + " | sed '1,/^-/d' > new.dat")

print("dc file extraction Done!")

os.system("awk '{print $28}' new.dat > total.dat")
os.system("awk '{print $2}' " + args.pdbfile + " > atomNo.dat")
os.system("awk '{print $6}' " + args.pdbfile + " > residue.dat")

print("file cutting Done!")

R_1 = open('total.dat',encoding='utf-8')
DcTotal = np.genfromtxt(R_1, dtype=np.float64, usecols=[0])
R_2 = open('residue.dat',encoding='utf-8')
ResidueNo = np.genfromtxt(R_2, dtype=np.float64, usecols=[0])
NoList = list(range(1, len(DcTotal) + 1))
DcList = []

for i in range(0, len(ResidueNo)):
    for j in range(0, len(NoList)):
        if int(ResidueNo[i]) == NoList[j]:
            if ResidueNo[i] % 1 == 0:
                DcList.append(round(DcTotal[j],2))
            if ResidueNo[i] % 1 != 0:
                DcList.append(0.00)
                
print("DC list generation Done!")
print("list comparing...")       
         
spaceNo = list()
residue = list()
count = 0
f = open("residue.dat", encoding='utf-8')
for line in f.readlines():
    count += 1
    residue.append(line)
    if line.strip() == '':
        spaceNo.append(count)
        
for i in range(0,len(spaceNo)):
    DcList.insert(spaceNo[i], 0.0)
    
while len(DcList) < len(residue):
    DcList.append(0.0)
    
print("list comparison Done!")

colspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26),
            (26, 27), (30, 38), (38, 46), (46, 54), (54, 60), (60, 66), (76, 78),
            (78, 80)]

names = ['ATOM', 'serial', 'name', 'altloc', 'resname', 'chainid', 'resseq',
         'icode', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'element', 'charge']

pdb = pd.read_fwf(args.pdbfile, names=names, colspecs=colspecs)
new_pdb = pdb.drop(columns=['tempfactor', 'element', 'charge'], axis=1)
new_pdb['tempfactor'] = DcList
new_pdb.fillna('', inplace=True)

atom = new_pdb['ATOM'].values.tolist()
atomNumber = new_pdb['serial'].values.tolist()
atomName = new_pdb['name'].values.tolist()
atomAltloc = new_pdb['altloc'].values.tolist()
residueName = new_pdb['resname'].values.tolist()
chainId = new_pdb['chainid'].values.tolist()
residueNum = new_pdb['resseq'].values.tolist()
icode = new_pdb['icode'].values.tolist()
x = new_pdb['x'].values.tolist()
y = new_pdb['y'].values.tolist()
z = new_pdb['z'].values.tolist()
occupancy = new_pdb['occupancy'].values.tolist()
tempfactor = new_pdb['tempfactor'].values.tolist()

print("PDB generation and cut Done!")

# "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(...)

print("saving as pdb...")

try:
    with open(args.output, 'w') as f:
        for i in range(0, len(atom)):
    #print("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format(atom[i], int(atomNumber[i]), atomName[i], atomAltloc[i], residueName[i], 
                                                    #chainId[i], int(residueNum[i]), icode[i], x[i], y[i], z[i], occupancy[i], tempfactor[i]))
            if occupancy[i] == '':
                pass
            else:
                f.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format(atom[i], int(atomNumber[i]), atomName[i], atomAltloc[i], residueName[i], 
                                                                                                        chainId[i], int(residueNum[i]), icode[i], x[i], y[i], z[i], occupancy[i], tempfactor[i]))
                f.write('\n')
except:
    print("error")
else:
    print("finished successfully")

os.system("rm total.dat")
os.system("rm residue.dat")
os.system("rm new.dat")
os.system("rm atomNo.dat")
