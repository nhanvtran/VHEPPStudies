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
files=glob.glob("/Users/ntran/Documents/Research/Ext/VHEPP/data/pgun_pi_1to10000gev%rfull005/pgun_pi1000gev*.slcio")
factory = LCFactory.getInstance()
reader = factory.createLCReader()

fileOutStatus1  = open('dat/of_status1.dat','w');
fileOutStatus3  = open('dat/of_status3.dat','w');
fileOutPanPFA   = open('dat/of_PanPFA.dat','w');
fileOutClusters = open('dat/of_CaloHits.dat','w');
fileOutTracks   = open('dat/of_Tracks.dat','w');
fileOutHCaloR = open('dat/of_HCaloHits_r.dat','w');
fileOutECaloR = open('dat/of_ECaloHits_r.dat','w');

fileOutSinglePion = open('dat/of_singlePi.dat','w');


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
		colCl = evt.getCollection("ReconClusters");
		colTr = evt.getCollection("Tracks");
		colHCB = evt.getCollection("HcalBarrelHits");
		colECB = evt.getCollection("EcalBarrelHits");

		nMc=col.getNumberOfElements()
		nPF=colPF.getNumberOfElements();
		nCl=colCl.getNumberOfElements();
		nTr=colTr.getNumberOfElements();
		nHCB=colHCB.getNumberOfElements();
		nECB=colECB.getNumberOfElements();		
		print "----------- ------------ ",nEvent
		# print "number of particles = ", nMc, nPF, nTr
		fileOutStatus1.write( '-99 0 0 0 0 \n' );
		fileOutStatus3.write( '-99 0 0 0 0 \n' );
		fileOutPanPFA.write( '-99 0 0 0 0 \n' );
		fileOutClusters.write( '-99 0 0 0 0 \n' );
		fileOutTracks.write( '-99 0 0 0 0 \n' );
		for i in range(nMc): # loop over all particles 
			par=col.getElementAt(i)
			if par.getGeneratorStatus() == 1: 

				singlePiOutput = '';
				maxsinglePiEnergy = -99;
				for i in range(nPF): # loop over all particles
					parPF=colPF.getElementAt(i) 
					# fileOutPanPFA.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} \n'.format(parPF.getType(),parPF.getMomentum()[0],parPF.getMomentum()[1],parPF.getMomentum()[2],parPF.getEnergy()) )
					trackp = -99;
					curstring = "energies [pfo,cluster,track]:"
					curstring += str(round(parPF.getEnergy(),4)) + " ";
					if parPF.getClusters().size() > 0: curstring += str(round(parPF.getClusters()[0].getEnergy(),4)) + " ";
					else: curstring += " - ";
					if parPF.getClusters().size() > 1: print "!!!!!!!!!!!!!!!!!!!!"
					if parPF.getTracks().size() > 0: 
							#http://flc.desy.de/lcnotes/notes/localfsExplorer_read?currentPath=/afs/desy.de/group/flc/lcnotes/LC-DET-2006-004.pdf
							trackp = math.fabs( (3e-4*5/parPF.getTracks()[0].getOmega())*math.sqrt(1 + parPF.getTracks()[0].getTanLambda()*parPF.getTracks()[0].getTanLambda() ) );
							curstring += str(round(math.fabs(trackp),4)) + " ";
					else: curstring += " - ";
					print curstring, "--- nclusters = ", str(parPF.getClusters().size()),",hits per cluster: ", parPF.getClusters()[0].getCalorimeterHits().size()

					if parPF.getEnergy() > maxsinglePiEnergy: 
						singlePiOutput = '{0:4.4f} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} {5:4.4f} {6:4.4f} {7:4.4f} {8:4.4f} {9:4.4f} {10:4.4f} \n'.format(par.getMomentum()[0],par.getMomentum()[1],par.getMomentum()[2],par.getEnergy(),parPF.getMomentum()[0],parPF.getMomentum()[1],parPF.getMomentum()[2],parPF.getEnergy(),parPF.getEnergy(), parPF.getClusters()[0].getEnergy(), trackp ) ;

						# ECluster = -99;
						# EClusterSummed = 0;
						# if parPF.getClusters().size() > 0: 
						# 	ECluster = parPF.getClusters()[0].getEnergy();
						# 	for i in range(parPF.getClusters()[0].getCalorimeterHits().size()):
						# 		EClusterSummed = EClusterSummed + parPF.getClusters()[0].getCalorimeterHits()[i].getEnergy()
						# 		#print parPF.getClusters()[0].getClusters().size(), i,":",parPF.getClusters()[0].getCalorimeterHits()[i].getType(),parPF.getClusters()[0].getCalorimeterHits()[i].getEnergy(),parPF.getClusters()[0].getCalorimeterHits()[i].getEnergyError();
						# #print "SummedEnergy = ", EClusterSummed
						
						# PTrack = -99;
						# if parPF.getTracks().size() > 0:
						# 	PTrack = (3e-4*5/parPF.getTracks()[0].getOmega())*math.sqrt(1 + parPF.getTracks()[0].getTanLambda()*parPF.getTracks()[0].getTanLambda() );
						# fileOutPanPFA.write( '{0:3} {1:4.4f} {2:4.4f} {3:4.4f} {4:4.4f} {5:4.4f} {6:4.4f} \n'.format(parPF.getType(),parPF.getMomentum()[0],parPF.getMomentum()[1],parPF.getMomentum()[2],parPF.getEnergy(),ECluster,PTrack ));

				fileOutSinglePion.write(singlePiOutput);



		if nEvent > 1000000: break;
reader.close() # close the file



