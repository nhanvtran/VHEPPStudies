import ROOT
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.10);
ROOT.gStyle.SetPalette(1);
import math

def makeCanvas(hists, tags):

	colors = [1,2,4,5,6,7]

	leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	tmax = -999;
	for i in range(len(hists)):
		leg.AddEntry(hists[i],tags[i],'l');
		if hists[i].GetMaximum() > tmax: tmax = hists[i].GetMaximum();

	c = ROOT.TCanvas("c","c",1000,800);
	hists[0].SetMinimum(0);
	hists[0].SetMaximum(tmax*1.2);
	for i in range(len(hists)):
		if i == 0: hists[i].Draw();
		else: 
			hists[i].SetLineColor( colors[i] );
			hists[i].Draw('sames');

	leg.Draw();
	c.SaveAs("plots/"+hists[0].GetName()+".pdf")
	c.SaveAs("plots/"+hists[0].GetName()+".png")

def getHists(tin,postfix):

	h_jpt      = ROOT.TH1F("h_jpt"+postfix,"; pT; N",100,0,10000);
	h_jp       = ROOT.TH1F("h_jp"+postfix,"; p; N",100,0,10000);
	h_jeta     = ROOT.TH1F("h_jeta"+postfix,"; eta; N",50,-5,5);
	h_jphi     = ROOT.TH1F("h_jphi"+postfix,"; phi; N",25,0,2*3.15);
	h_jmass    = ROOT.TH1F("h_jmass"+postfix,"; mass; N",50,0,200);
	h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,"; soft drop mass (#beta = 0); N",50,0,200);
	h_njets    = ROOT.TH1F("h_njets"+postfix,"; njets; N",10,0,10);

	# h_jpt      = ROOT.TH1F("h_jpt"+postfix,"; pT; N",100,0,1000);
	# h_jp       = ROOT.TH1F("h_jp"+postfix,"; p; N",100,0,1000);
	# h_jeta     = ROOT.TH1F("h_jeta"+postfix,"; eta; N",50,-5,5);
	# h_jphi     = ROOT.TH1F("h_jphi"+postfix,"; phi; N",25,0,2*3.15);
	# h_jmass    = ROOT.TH1F("h_jmass"+postfix,"; mass; N",50,0,200);
	# h_jmass_sd = ROOT.TH1F("h_jmass_sd"+postfix,"; soft drop mass (#beta = 0); N",50,0,200);
	# h_njets    = ROOT.TH1F("h_njets"+postfix,"; njets; N",10,0,10);

	# print tin.GetEntriesFast();
	for i in range(tin.GetEntriesFast()):
		tin.GetEntry(i);
		for j in range(len(tin.jpt)):
			# print tin.jisleptag[j]
			if tin.jisleptag[j] == 0 and math.fabs(tin.jeta[j]) < 1.0: 
				# print tin.jpt[j]
				h_jpt.Fill( tin.jpt[j] );
				h_jp.Fill( tin.jp[j] );				
				h_jeta.Fill( tin.jeta[j] );
				h_jphi.Fill( tin.jphi[j] );
				h_jmass.Fill( tin.jmass[j] );
				h_jmass_sd.Fill( tin.jmass_sd[j] );	
		h_njets.Fill( tin.njets );

	hists = [];
	hists.append( h_jpt );
	hists.append( h_jp );
	hists.append( h_jeta );
	hists.append( h_jphi );
	hists.append( h_jmass );
	hists.append( h_jmass_sd );
	hists.append( h_njets );

	return hists;

def MCinfo(tin):

	h_mzp = ROOT.TH1F("h_mzp","; zprime mass; N",50,9000,11000);
	for i in range(tin.GetEntriesFast()):
		tin.GetEntry(i);
		h_mzp.Fill(tin.gen_mZp);
	return h_mzp;

def rms90(h):

	nbins = h.GetXaxis().GetNbins();
	imean = h.FindBin( h.GetMean() );
	entries = 0.9*h.GetEntries();
	w = h.GetBinContent( imean );
	x = h.GetBinCenter( imean );
	sumw = w;
	sumwx = w*x;
	sumwx2 = w*x*x;
	for i in range(1,nbins):
		if (i>0):
			w = h.GetBinContent( imean - i );
			x = h.GetBinCenter( imean - i );
			sumw += w;
			sumwx += w*x;
			sumwx2 += w*x*x;
		if (i<=nbins):
			w = h.GetBinContent( imean - i );
			x = h.GetBinCenter( imean - i );
			sumw += w;
			sumwx += w*x;
			sumwx2 += w*x*x;
		if sumw > entries: break;

	x = sumwx/sumw;
	rms90val = math.fabs(sumwx2/sumw - x*x);
	print "rms90 = ", rms90val, h.GetRMS();
	return rms90;

def GetResolutions(tin1,tin2,pf):

	h_jpres = ROOT.TH1F("h_jpres_"+pf, "; (P - PGEN)/PGEN; au", 50,-1.0,1.0);

	jpres_vals = [];
	for i in range(tin1.GetEntriesFast()):
		tin1.GetEntry(i);
		tin2.GetEntry(i);

		for j in range(len(tin1.jpt)):
			if tin1.jisleptag[j] == 0 and math.fabs(tin1.jeta[j]) < 0.4: 
				for k in range(len(tin2.jpt)):

					dr = math.sqrt( (tin1.jphi[j]-tin2.jphi[k])*(tin1.jphi[j]-tin2.jphi[k]) + (tin1.jeta[j]-tin2.jeta[k])*(tin1.jeta[j]-tin2.jeta[k]) );
					if dr < 0.4: 
						h_jpres.Fill( (tin1.jp[j] - tin2.jp[k])/tin2.jp[k] );		
						jpres_vals.append( (tin1.jp[j] - tin2.jp[k])/tin2.jp[k] );
	
	jpres_vals.sort();
	loit = int(0.05*len(jpres_vals));
	hiit = int(0.95*len(jpres_vals));
	h_jpres_rms90 = ROOT.TH1F("h_jpres_rms90"+pf, "; (P - PGEN)/PGEN; au", 500,-1,1);
	for v in range(loit,hiit):
		h_jpres_rms90.Fill( jpres_vals[v] );
	rms90val = h_jpres_rms90.GetRMS();
	return (h_jpres,rms90val);

if __name__ == '__main__':

	fa = ROOT.TFile("dat/of_PanPFA.root");
	ta = fa.Get("tPFA");
	tb = fa.Get("tGEN");
	tc = fa.Get("tcalo");
	tmc = fa.Get("tMC");


	print "Getting hists..."
	ha = getHists(ta,'a');
	hb = getHists(tb,'b');
	hc = getHists(tc,'c');	

	for i in range(len(ha)):
		makeCanvas( [ha[i],hb[i],hc[i]], ['pf','gen','calo'] )

	hmc = MCinfo(tmc);
	makeCanvas([hmc],['particlelevel'])
	# hres = Resolutions(tb,ta,tc);

	hres_pf,rms90_pf = GetResolutions( ta, tb, "pf" );
	hres_cal,rms90_cal = GetResolutions( tc, tb, "calo" );
	makeCanvas( [hres_pf,hres_cal], ['pf,rms90 ='+str(round(rms90_pf,4)),'calo,rms90 ='+str(round(rms90_cal,4))] )




