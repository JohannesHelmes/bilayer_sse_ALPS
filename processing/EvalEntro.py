import matplotlib.pylab as plt
import pyalps,argparse
import os
import plottery
import uncertainties as unc
import numpy as np

parser = argparse.ArgumentParser(description='Evaluate Renyi entropies for range of L', epilog='(C) Johannes Helmes 2014')

parser.add_argument('--infile','-i', help='Prefix of result files',required=True)
parser.add_argument('--IncNoRange','-n', nargs=2,type=int,help='range of IncNo steps: start stop (int, int)')
parser.add_argument('--foreach','-f',default='h',help='Parameter name, (default h)')
parser.add_argument('--plot','-p',action='store_true')
parser.add_argument('--verbose','-v',action='store_true')
parser.add_argument('--partitioned','-a',help="set this flag, if partitioned run",action='store_true')
args=parser.parse_args()

REntropy={}
if (args.IncNoRange!=None):
    IncNos=range(args.IncNoRange[0],args.IncNoRange[1]+1)
    #IncNos = [%.1f % i for i in range(args.IncNoRange[0], args.IncNoRange[1])]
else:
    IncNos=range(100)

data=[]
if args.partitioned:
    i=0
    while (os.path.isfile(args.infile+"."+str(i)+".out.xml")):
        data = data+pyalps.loadMeasurements(pyalps.getResultFiles(prefix=args.infile+"."+str(i)),['ED','EG'])
        i+=1
else:
    data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix=args.infile),['ED','EG'])


if args.verbose:
    print IncNos
    print data

#renyi_dataD = pyalps.collectXY(data, x='IncNo', y='ED', foreach=['h'])
renyi_dataG = pyalps.collectXY(data, x='IncNo', y='EG', foreach=[args.foreach])


if args.verbose:
    print renyi_dataG


result_set_x,result_set_y,result_set_yerr=plottery.renyi_sse_add(renyi_dataG, for_each=str(args.foreach),inc_name='IncNo',inc_range=IncNos)
plottery.print_data(result_set_x,result_set_y,result_set_yerr)

if args.plot:
#plt.ylim([0, 2])
    plottery.empty_kiln()
    plt.ylabel("$S_2(A)$")
    plt.xlabel(args.foreach)
    plottery.errorbar(result_set_x,result_set_y,yerr=result_set_yerr, fmt='o')
    plottery.fire(type='sequential')
    plt.show()
#plt.savefig("/home/helmes/Desktop/test_johanes.pdf")
