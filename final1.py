"""
Back up any previous excel spreadsheets. Running this simulation
Will overwrite any existing data
Written to file
"""

from prody import *
import pylab as op
op.ion()
protein=parsePDB('2DX3')
import random
import time
import numpy
import matplotlib
import xlsxwriter as excel
import math
#instant.
indmolelist=[(1.448680*(10**(-15))),(1.448680*(10**(-18))),(1.448680*(10**(-21))),(1.448680*(10**(-24))),(1.448680*(10**(-27)))] 
avagadro=6.022*(10**23)


#prelim excel writing
#spreadsheet=excel.Workbook('Data.xlsx')
#rawdatasheet=spreadsheet.add_worksheet('rawdata')
#rawdatasheet.write(0, 0, "Peptide Atom")
#rawdatasheet.write(0, 1, "Searched Atom")
#rawdatasheet.write(0,2,"x")
#rawdatasheet.write(0,3,"y")
#rawdatasheet.write(0,4,"z")
#rawdatasheet.write(0,5,"searched residue")
#rawdatasheet.write(0,6,"searching residue")
excel.Workbook('1.csv').close()
excel.Workbook('2.csv').close()
excel.Workbook('3.csv').close()
excel.Workbook('4.csv').close()
excel.Workbook('5.csv').close()



#basic PDB parse
sequence=protein.getResnames()
coordinates=protein.getCoords()
coordinates=coordinates.tolist()
atoms=Atom.copy(protein)
absoluteAtomNames=[]

#getting absolute atom names
for all2 in range (0,len(atoms)):
    current=str(atoms[all2])
    current=current.split(" ")
    toappend=current[1][0]
    absoluteAtomNames.append(toappend)

#appending to single lists to calculate poles
x=[]
y=[]
z=[]
for p in range (0,len(coordinates)):
    x.append(coordinates[p][0])
    y.append(coordinates[p][1])
    z.append(coordinates[p][2])

#finding min at element zero via sort
x.sort()
y.sort()
z.sort()

#making simulation environnment/box
xmin=x[0]-4
xmax=x[-1]+4
ymin=y[0]-4
ymax=y[-1]+4
zmin=z[0]-4
zmax=z[-1]+4
print (xmin, xmax)

#establishing variables 
totalvolume=(abs(xmin-xmax)*abs(ymin-ymax)*abs(zmin-zmax)) #total volume
specificvolume=(4**3) #box volume for op.
ratio=(specificvolume/totalvolume)

print (ratio)
print (totalvolume)
print (specificvolume)
print (abs(xmin-xmax))
print (abs(ymin-ymax))
print (abs(zmin-zmax))

#iterations and writes
for simulationnumber in range(0,5):
    fileString=str(simulationnumber+1)+".csv"
    datafile=open(fileString,'w')
    datafile.close()
    xyzlist=[]
    atombindinglist=[]
    atomsearchinglist=[]
    residuebindinglist=[]
    residuesearchinglist=[]
    distlist=[]
    totalsimulationmolecules=(avagadro*(indmolelist[simulationnumber]))
    moleculesinspecificarea=(ratio*totalsimulationmolecules)
    oxygenorientationcof=0.5
    hydrogenorientationcof=0.15937638889
    for peptidemolecule in range (0,len(absoluteAtomNames)):
        print (peptidemolecule)
        searchingresidue=sequence[peptidemolecule]
        atom=absoluteAtomNames[peptidemolecule]
        if atom == "N":
            quantityiterations=(int(moleculesinspecificarea*hydrogenorientationcof))
            for indiceneverused in range (0,quantityiterations):
                xyzcurrent=(random.uniform(float(coordinates[peptidemolecule][0])-2,float(coordinates[peptidemolecule][0])+2),random.uniform(float(coordinates[peptidemolecule][1])-2, float(coordinates[peptidemolecule][1])+2),random.uniform(float(coordinates[peptidemolecule][2])-2,float(coordinates[peptidemolecule][2])+2))
                bindatom=("H")
                searchingresidue=("H2O")
                xyzlist.append(xyzcurrent)
                atombindinglist.append(bindatom)
                residuebindinglist.append(searchingresidue)
                residuesearchinglist.append(sequence[peptidemolecule])
                atomsearchinglist.append(atom)
            
            for atompeptideindice in range (0,len(absoluteAtomNames)):
                if absoluteAtomNames[atompeptideindice] == "H" and sequence[atompeptideindice] != sequence[peptidemolecule]:
                    xyzcurrent=(coordinates[atompeptideindice][0],coordinates[atompeptideindice][1],coordinates[atompeptideindice][2])
                    bindatom=("H")
                    searchingresidue=sequence[atompeptideindice]
                    xyzlist.append(xyzcurrent)
                    atombindinglist.append(bindatom)
                    residuebindinglist.append(searchingresidue)
                    residuesearchinglist.append(sequence[peptidemolecule])
                    atomsearchinglist.append(atom)
            datafile=open(fileString,'a')
            iterator=0
            for iterator in range (0,len(xyzlist)-1):
                x=(xyzlist[iterator][0])
                y=(xyzlist[iterator][1])
                z=(xyzlist[iterator][2])
                distance=math.sqrt((abs(x-coordinates[peptidemolecule][0]))**2+((abs(y-coordinates[peptidemolecule][1])))**2+((abs(z-coordinates[peptidemolecule][2])**2)))
                distlist.append(distance)
            distlist,xyzlist,atombindinglist,atomsearchinglist,residuebindinglist,residuesearchinglist=zip(*sorted(zip(distlist,xyzlist,atombindinglist,atomsearchinglist,residuebindinglist,residuesearchinglist)))


            x=(str(xyzlist[0][0]))
            y=(str(xyzlist[0][1]))
            z=(str(xyzlist[0][2]))
            o=(str(atombindinglist[0]))
            p=(str(atomsearchinglist[0]))
            w=(str(residuebindinglist[0]))
            i=(str(residuesearchinglist[0]))
            j=(str(coordinates[peptidemolecule][0]))
            k=(str(coordinates[peptidemolecule][1]))
            l=(str(coordinates[peptidemolecule][2]))
            toappend=(x+","+y+","+z+","+ o+","+p+","+w+","+i+","+j+","+k+","+l+"\n")
            datafile.write(toappend)


            xyzlist=[]
            atombindinglist=[]
            atomsearchinglist=[]
            residuebindinglist=[]
            residuesearchinglist=[]
            distlist=[]
            datafile.close()

        elif atom == "C":
            quantityiterations=(int(moleculesinspecificarea*hydrogenorientationcof))  #also favourable hydrogen orientations
            for indiceneverused in range (0,quantityiterations):
                xyzcurrent=(random.uniform(float(coordinates[peptidemolecule][0])-2,float(coordinates[peptidemolecule][0])+2),random.uniform(float(coordinates[peptidemolecule][1])-2, float(coordinates[peptidemolecule][1])+2),random.uniform(float(coordinates[peptidemolecule][2])-2,float(coordinates[peptidemolecule][2])+2))
                bindatom=("H")
                searchingresidue=("H2O")
                xyzlist.append(xyzcurrent)
                atombindinglist.append(bindatom)
                residuebindinglist.append(searchingresidue)
                residuesearchinglist.append(sequence[peptidemolecule])
                atomsearchinglist.append(atom)
            
            for atompeptideindice in range (0,len(absoluteAtomNames)):
                if absoluteAtomNames[atompeptideindice] == "H" and sequence[atompeptideindice] != sequence[peptidemolecule]:
                    xyzcurrent=(coordinates[atompeptideindice][0],coordinates[atompeptideindice][1],coordinates[atompeptideindice][2])
                    bindatom=("H")
                    searchingresidue=sequence[atompeptideindice]
                    xyzlist.append(xyzcurrent)
                    atombindinglist.append(bindatom)
                    residuebindinglist.append(searchingresidue)
                    residuesearchinglist.append(sequence[peptidemolecule])
                    atomsearchinglist.append(atom)

            datafile=open(fileString,'a')
            for iterator in range (0,len(xyzlist)):
                x=(xyzlist[iterator][0])
                y=(xyzlist[iterator][1])
                z=(xyzlist[iterator][2])
                distance=math.sqrt((abs(x-coordinates[peptidemolecule][0]))**2+((abs(y-coordinates[peptidemolecule][1])))**2+((abs(z-coordinates[peptidemolecule][2])**2)))
                distlist.append(distance)
            distlist,xyzlist,atombindinglist,atomsearchinglist,residuebindinglist,residuesearchinglist=zip(*sorted(zip(distlist,xyzlist,atombindinglist,atomsearchinglist,residuebindinglist,residuesearchinglist)))


            x=(str(xyzlist[0][0]))
            y=(str(xyzlist[0][1]))
            z=(str(xyzlist[0][2]))
            o=(str(atombindinglist[0]))
            p=(str(atomsearchinglist[0]))
            w=(str(residuebindinglist[0]))
            i=(str(residuesearchinglist[0]))
            j=(str(coordinates[peptidemolecule][0]))
            k=(str(coordinates[peptidemolecule][1]))
            l=(str(coordinates[peptidemolecule][2]))
            toappend=(x+","+y+","+z+","+ o+","+p+","+w+","+i+","+j+","+k+","+l+"\n")
            datafile.write(toappend)


            xyzlist=[]
            atombindinglist=[]
            atomsearchinglist=[]
            residuebindinglist=[]
            residuesearchinglist=[]
            distlist=[]
            datafile.close()
        
        
        elif atom == "O":
            quantityiterations=(int(moleculesinspecificarea*hydrogenorientationcof))  #also favourable hydrogen orientations
            for indiceneverused in range (0,quantityiterations):
                xyzcurrent=(random.uniform(float(coordinates[peptidemolecule][0])-2,float(coordinates[peptidemolecule][0])+2),random.uniform(float(coordinates[peptidemolecule][1])-2, float(coordinates[peptidemolecule][1])+2),random.uniform(float(coordinates[peptidemolecule][2])-2,float(coordinates[peptidemolecule][2])+2))
                bindatom=("H")
                searchingresidue=("H2O")
                xyzlist.append(xyzcurrent)
                atombindinglist.append(bindatom)
                residuebindinglist.append(searchingresidue)
                residuesearchinglist.append(sequence[peptidemolecule])
                atomsearchinglist.append(atom)
            
            for atompeptideindice in range (0,len(absoluteAtomNames)):
                if absoluteAtomNames[atompeptideindice] == "H" and sequence[atompeptideindice] != sequence[peptidemolecule]:
                    xyzcurrent=(coordinates[atompeptideindice][0],coordinates[atompeptideindice][1],coordinates[atompeptideindice][2])
                    bindatom=("H")
                    searchingresidue=sequence[atompeptideindice]
                    xyzlist.append(xyzcurrent)
                    atombindinglist.append(bindatom)
                    residuebindinglist.append(searchingresidue)
                    residuesearchinglist.append(sequence[peptidemolecule])
                    atomsearchinglist.append(atom)

            datafile=open(fileString,'a')
            for iterator in range (0,len(xyzlist)):
                x=(xyzlist[iterator][0])
                y=(xyzlist[iterator][1])
                z=(xyzlist[iterator][2])
                distance=math.sqrt((abs(x-coordinates[peptidemolecule][0]))**2+((abs(y-coordinates[peptidemolecule][1])))**2+((abs(z-coordinates[peptidemolecule][2])**2)))
                distlist.append(distance)
            distlist,xyzlist,atombindinglist,atomsearchinglist,residuebindinglist,residuesearchinglist=zip(*sorted(zip(distlist,xyzlist,atombindinglist,atomsearchinglist,residuebindinglist,residuesearchinglist)))


            x=(str(xyzlist[0][0]))
            y=(str(xyzlist[0][1]))
            z=(str(xyzlist[0][2]))
            o=(str(atombindinglist[0]))
            p=(str(atomsearchinglist[0]))
            w=(str(residuebindinglist[0]))
            i=(str(residuesearchinglist[0]))
            j=(str(coordinates[peptidemolecule][0]))
            k=(str(coordinates[peptidemolecule][1]))
            l=(str(coordinates[peptidemolecule][2]))
            toappend=(x+","+y+","+z+","+ o+","+p+","+w+","+i+","+j+","+k+","+l+"\n")
            datafile.write(toappend)


            xyzlist=[]
            atombindinglist=[]
            atomsearchinglist=[]
            residuebindinglist=[]
            residuesearchinglist=[]
            distlist=[]
            datafile.close()
        elif atom == "H":
            quantityiterations=(int(moleculesinspecificarea*0.5))
            for indiceneverused in range (0,quantityiterations):
                xyzcurrent=(random.uniform(float(coordinates[peptidemolecule][0])-2,float(coordinates[peptidemolecule][0])+2),random.uniform(float(coordinates[peptidemolecule][1])-2, float(coordinates[peptidemolecule][1])+2),random.uniform(float(coordinates[peptidemolecule][2])-2,float(coordinates[peptidemolecule][2])+2))
                bindatom=("O")
                searchingresidue=("H2O")
                xyzlist.append(xyzcurrent)
                atombindinglist.append(bindatom)
                residuebindinglist.append(searchingresidue)
                residuesearchinglist.append(sequence[peptidemolecule])
                atomsearchinglist.append(atom)
            
            for atompeptideindice in range (0,len(absoluteAtomNames)):
                if absoluteAtomNames[atompeptideindice] == "O" and sequence[atompeptideindice] != sequence[peptidemolecule]:
                    xyzcurrent=(coordinates[atompeptideindice][0],coordinates[atompeptideindice][1],coordinates[atompeptideindice][2])
                    bindatom=("O")
                    searchingresidue=sequence[atompeptideindice]
                    xyzlist.append(xyzcurrent)
                    atombindinglist.append(bindatom)
                    residuebindinglist.append(searchingresidue)
                    residuesearchinglist.append(sequence[peptidemolecule])
                    atomsearchinglist.append(atom)
               
                elif absoluteAtomNames[atompeptideindice] == "N" and sequence[atompeptideindice] != sequence[peptidemolecule]:
                    xyzcurrent=(coordinates[atompeptideindice][0],coordinates[atompeptideindice][1],coordinates[atompeptideindice][2])
                    bindatom=("N")
                    searchingresidue=sequence[atompeptideindice]
                    xyzlist.append(xyzcurrent)
                    atombindinglist.append(bindatom)
                    residuebindinglist.append(searchingresidue)
                    residuesearchinglist.append(sequence[peptidemolecule])
                    atomsearchinglist.append(atom)
           
            datafile=open(fileString,'a')
            for iterator in range (0,len(xyzlist)):
                x=(xyzlist[iterator][0])
                y=(xyzlist[iterator][1])
                z=(xyzlist[iterator][2])
                distance=math.sqrt((abs(x-coordinates[peptidemolecule][0]))**2+((abs(y-coordinates[peptidemolecule][1])))**2+((abs(z-coordinates[peptidemolecule][2])**2)))
                distlist.append(distance)
            distlist,xyzlist,atombindinglist,atomsearchinglist,residuebindinglist,residuesearchinglist=zip(*sorted(zip(distlist,xyzlist,atombindinglist,atomsearchinglist,residuebindinglist,residuesearchinglist)))


            x=(str(xyzlist[0][0]))
            y=(str(xyzlist[0][1]))
            z=(str(xyzlist[0][2]))
            o=(str(atombindinglist[0]))
            p=(str(atomsearchinglist[0]))
            w=(str(residuebindinglist[0]))
            i=(str(residuesearchinglist[0]))
            j=(str(coordinates[peptidemolecule][0]))
            k=(str(coordinates[peptidemolecule][1]))
            l=(str(coordinates[peptidemolecule][2]))
            toappend=(x+","+y+","+z+","+ o+","+p+","+w+","+i+","+j+","+k+","+l+"\n")
            datafile.write(toappend)


            xyzlist=[]
            atombindinglist=[]
            atomsearchinglist=[]
            residuebindinglist=[]
            residuesearchinglist=[]
            distlist=[]
            datafile.close()
            

        

    #writing everything to raw file and creating a new sheet

print (xyzlist)
print (atombindinglist)
print (residuebindinglist)
print (residuesearchinglist)
print (atomsearchinglist)
