# Example: read truth-level particles and create Z peak  
# Use http://atlaswww.hep.anl.gov/hepsim/show.php?item=146&reco=rfull001 
#  Get 4 files in 2 threads (see HepSim manual):
#  > hs-get gev250ee_pythia6_zpole_ee%rfull001 gev250ee_pythia6_zpole_ee 2 4 
# S.Chekanov (ANL)

from org.lcsim.lcio import LCIOReader
from hep.io.sio import SIOReader
from hep.lcio.implementation.sio import SIOLCReader
from hep.lcio.implementation.io import LCFactory
from hep.lcio.event import * 
from hep.lcio.io import *
from jhplot import *  # import graphics
from hephysics.particle import LParticle
import math

# make list of files..
import glob
files=glob.glob("dat/*.slcio")
factory = LCFactory.getInstance()
reader = factory.createLCReader()

fileOutStatus1 = open('dat/of_status1.dat','w');
fileOutStatus3 = open('dat/of_status3.dat','w');
fileOutPanPFA  = open('dat/of_PanPFA.dat','w');

nEvent=0
for f in files:
	print "Open file=",f
	reader.open(f) 
	while(1):
		evt=reader.readNextEvent()
		if (evt == None): break
		nEvent=nEvent+1
		# print " file event: ",evt.getEventNumber(), " run=",evt.getRunNumber()
		if (nEvent%100==0): print "# Event: ",nEvent
		col = evt.getCollection("MCParticle")
		colPF = evt.getCollection("PandoraPFOCollection");
		colTr = evt.getCollection("Tracks");
		nMc=col.getNumberOfElements()
		nPF=colPF.getNumberOfElements();
		nTr=colTr.getNumberOfElements();
		print " ----------- ------------ ",nEvent
		print "number of particles = ", nMc, nPF, nTr
		fileOutStatus1.write( '0 0 0 0 0 \n' );
		fileOutStatus3.write( '0 0 0 0 0 \n' );
		fileOutPanPFA.write( '0 0 0 0 0 \n' );
		for i in range(nMc): # loop over all particles 
			par=col.getElementAt(i)
			if par.getGeneratorStatus() == 3: 
				pdg=par.getPDG();
				charge=par.getCharge();
				if par.getParents().size() > 0:
					if math.fabs(par.getParents()[0].getPDG()) == 24:
						print i, "pdg=",pdg, ", status=",par.getGeneratorStatus(), " mass=",par.getMass(), " mother=",par.getParents()[0].getPDG(), " momsize=",len(par.getMomentum()); 
						fileOutStatus3.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(par.getPDG(),par.getMomentum()[0],par.getMomentum()[1],par.getMomentum()[2],par.getEnergy()) )
			if par.getGeneratorStatus() == 1: 
				if (math.fabs(par.getPDG()) == 12 or math.fabs(par.getPDG()) == 14 or math.fabs(par.getPDG()) == 16): continue;
				fileOutStatus1.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(par.getPDG(),par.getMomentum()[0],par.getMomentum()[1],par.getMomentum()[2],par.getEnergy()) )
		
		for i in range(nPF): # loop over all particles
			parPF=colPF.getElementAt(i) 
			fileOutPanPFA.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parPF.getType(),parPF.getMomentum()[0],parPF.getMomentum()[1],parPF.getMomentum()[2],parPF.getEnergy()) )


		if nEvent > 1000: break;

		# particles=[]
		# pp=0
		# nn=0
		# for i in range(nMc): # loop over all particles 
		# 	par=col.getElementAt(i)
		# 	if(par.getGeneratorStatus() == 1 and par.getCharge() !=0): 
		# 		vertex = par.getVertex();
		# 		pdg=par.getPDG()
		# 		momentum = par.getMomentum() 
		# 		ee=par.getEnergy()
		# 		mass=par.getMass()
		# 		px=momentum[0]
		# 		py=momentum[1]
		# 		pz=momentum[2]
		# 		charge=par.getCharge()
		# 		p=LParticle(px,py,pz,ee,mass) 
		# 		p.setCharge(charge)
		# 		particles.append(p)
		# 		print nn, "pdg=",pdg, " ",px,py,pz, " charge=",charge 
		# 		nn=nn+1

reader.close() # close the file



