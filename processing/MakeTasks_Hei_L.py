#!/bin/bash alpspython

import pyalps,sys,argparse,os
import numpy as np

parser = argparse.ArgumentParser(description='Create input files for ALPS simulations of Heisenberg models with constant L', epilog='(C) Johannes Helmes 2013')


parser.add_argument('--sweeps','-s', type=int,default=1000000,help='Number of sweeps, default=1000000')
parser.add_argument('--therm','-t', type=int,default=10000,help='Number of thermalization steps, default=10000')
parser.add_argument('--length','-l', type=int,help='Lattice size',required=True)
parser.add_argument('--beta','-b', type=float,help='beta-1/T, default=l')
parser.add_argument('--Jrange','-j', nargs=3,type=float,help='range of J values: start stop inc (float, float, float)')
parser.add_argument('--infile','-i', help='Prefix of .in.xml')
parser.add_argument('--geofile','-g', help='File (prefix) containing the geometry of the increments')
parser.add_argument('--partsize','-p',default=1, type=int,help='Size subsimulation partition for better distribution on HPC, default=1')
parser.add_argument('--half', help='Indicates, that only half of the geometry file is used (e.g. only up to the square)',action='store_true')
args=parser.parse_args()

if (args.Jrange==None):
    #Jlist=[0.0,0.5,1.0,1.5]+np.arange(2.0,3.6,0.1).tolist()
    #Jlist=[0.0,0.25,1.0,1.5,2.522,3.0]
    #Jlist=[5.46]
    #Jlist=[1.0,2.0,2.5,2.6,2.7,2.7755,2.8,2.9,3.0]
    Jlist=[0.0,1.0,2.0,3.0,4.0,5.0,5.1,5.2,5.3,5.4,5.46,5.5,5.6,5.7,5.8,5.9,6.0]
else:
    Jlist= np.arange(args.Jrange[0],args.Jrange[1],args.Jrange[2])
parms = []

beta=args.beta
if (beta==None):
    beta=args.length

for J in Jlist:
    Geometry=np.genfromtxt(args.geofile+'.geo',dtype=str)
    if (args.half):
        Geometry = Geometry[:len(Geometry)/2]

    for i,IncEl in enumerate(Geometry):
 
        parms.append(
            {
                'LATTICE'	: "bilayer",
                'LATTICE_LIBRARY' : "bilayerlib.xml",
                'Interaction' : 'A',
                'beta'		: beta,
                'THERMALIZATION': args.therm,
                'SWEEPS'	: args.sweeps,
                'Interlayer': J,
                'L'		: args.length,
                'M'     : 100,
                'IncStep': IncEl,
                'IncNo' : i
            }
        )

cuts=np.append(np.arange(0,len(parms),args.partsize),len(parms))
print cuts
path=os.path.dirname(args.infile)
for i,[s,e] in enumerate(zip(cuts[:-1],cuts[1:])):
    fname=path+'/'+os.path.basename(args.infile)+'.'+str(i)
    input_file = pyalps.writeInputFiles(fname,parms[s:e])
