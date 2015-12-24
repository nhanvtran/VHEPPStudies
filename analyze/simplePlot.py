import ROOT
import tdrstyle
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPadRightMargin(0.15);
ROOT.gStyle.SetPalette(1);


def makeCanvas(hists, tags):

	colors = [1,2,4,5,6,7]

	leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	for i in range(len(hists)):
		leg.AddEntry(hists[i],tags[i],'l');

	c = ROOT.TCanvas("c","c",1000,800);
	hists[0].SetMinimum(0);
	for i in range(len(hists)):
		if i == 0: hists[i].Draw();
		else: 
			hists[i].SetLineColor( colors[i] );
			hists[i].Draw('sames');

	leg.Draw();
	c.SaveAs("plots/"+hists[0].GetName()+".pdf")

def getHists(tin):

	h_jpt      = ROOT.TH1F("h_jpt","; pT; N",70,0,7000);
	h_jeta     = ROOT.TH1F("h_jeta","; eta; N",50,-5,5);
	h_jphi     = ROOT.TH1F("h_jphi","; phi; N",50,0,2*3.15);
	h_jmass    = ROOT.TH1F("h_jmass","; mass; N",50,0,200);
	h_jmass_sd = ROOT.TH1F("h_jmass_sd","; soft drop mass (#beta = 0); N",50,0,200);
	h_jpt_calo = ROOT.TH1F("h_jpt_calo","; pT; N",70,0,7000);

	print tin.GetEntriesFast();
	for i in range(tin.GetEntriesFast()):
		tin.GetEntry(i);
		for j in range(len(tin.jpt)):
			# print tin.jisleptag[j]
			if tin.jisleptag[j] == 0: 
				print tin.jpt[j]
				h_jpt.Fill( tin.jpt[j] );
				h_jeta.Fill( tin.jeta[j] );
				h_jphi.Fill( tin.jphi[j] );
				h_jmass.Fill( tin.jmass[j] );
				h_jmass_sd.Fill( tin.jmass_sd[j] );
		for j in range(len(tin.jpt_calo)):
			h_jpt_calo.Fill( tin.jpt_calo[j] );

	hists = [];
	hists.append( h_jpt );
	hists.append( h_jeta );
	hists.append( h_jphi );
	hists.append( h_jmass );
	hists.append( h_jmass_sd );
	hists.append( h_jpt_calo );
	return hists;

if __name__ == '__main__':

	fa = ROOT.TFile("dat/of_status1.root");
	fb = ROOT.TFile("dat/of_PanPFA.root");
	ta = fa.Get("t");
	tb = fb.Get("t");

	print "Getting hists..."
	ha = getHists(ta);
	hb = getHists(tb);

	for i in range(len(ha)):
		makeCanvas( [ha[i],hb[i]], ['gen','pf'] )

	makeCanvas( [ha[0],hb[0],ha[5]], ['gen','pf','calo'])


