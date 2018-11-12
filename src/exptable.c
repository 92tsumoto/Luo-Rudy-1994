#include "syspara.h"

void make_ExpTable()
{

	int vindex,kiindex;
	double v,ki;
	double am,bm,ah,bh,aj,bj;
	double ar,br,at,bt;
   	double ad,bd,af,bf;
   	double ax,bx;
   	double axs,bxs;

	for(vindex=0;vindex<VNMAX;vindex++){

        //v = (double)vindex/dvm-200.0;
		v = (double)vindex/dvm-(double)Emax;
		
		if(fabs(v+47.13)>0.001){
			am = 0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
		} else {
			am = 3.2;
		}
		bm = 0.08*exp(-v/11.0);

		ina.Tmss[vindex] = am/(am+bm);
		ina.Ttaum[vindex] = 1.0/(am+bm);

		if(v < -40.0){
			ah = 0.135*exp((80.0+v)/-6.8);
			bh = 3.56*exp(0.079*v)+3.1E+5*exp(0.35*v);
			aj = (-127140.0*exp(0.244*v)-3.474E-5*exp(-0.04391*v))*(v+37.78)/(1.0+exp(0.311*(v+79.23)));
			bj = (0.1212*exp(-0.01052*v))/(1.0+exp(-0.1378*(v+40.14)));
		} else {
			ah = 0.0;
			bh = 1.0/(0.13*(1.0+exp(-(v+10.66)/11.1)));
			aj = 0.0;
			bj = (0.3*exp(-2.535E-7*v))/(1.0+exp(-0.1*(v+32.0)));
		}
		ina.Thss[vindex] = ah/(ah+bh);
		ina.Ttauh[vindex] = 1.0/(ah+bh);
		ina.Tjss[vindex] = aj/(aj+bj);
		ina.Ttauj[vindex] = 1.0/(aj+bj);

		// for ical
		ad = 1.0/(1.0+exp(-(v+10.0)/6.24)); 
		if(fabs(v+10.0)>0.001){
			bd = ad*(1.0-exp(-(v+10.0)/6.24))/(0.035*(v+10.0));
		} else {
			bd = ad*(6.24*0.035);
		}
		ical.Tdss[vindex] = ad;
		ical.Ttaud[vindex] = bd;
		
		ical.Tfss[vindex] = 1.0/(1.0+exp((v+32.0)/8.0)) + 0.6/(1.0+exp((50.0-v)/20.0));
		ical.Ttauf[vindex] = 1.0/(0.0197*exp(-(0.0337*(v+10.0))*(0.0337*(v+10.0)))+0.02);

		// for ikx
		if(fabs(v+30.0)>0.001){
			ax = 7.19E-5*(v+30.0)/(1.0-exp(-0.148*(v+30.0)));
		}else {
			ax = 7.19E-5/0.148; 
		}
		if(fabs(v+30.0)>0.001){
			bx = 1.31E-4*(v+30.0)/(-1.0+exp(0.0687*(v+30.0)));
		} else {
			bx = 1.31E-4/30.0;
		}

		ikx.Txss[vindex] = ax/(ax+bx);
		ikx.Ttaux[vindex] = 1.0/(ax+bx);
		ikx.Trik[vindex] = 1.0/(1.0+exp((v-56.26)/32.1));

		// ikp
		ikp.Txkp[vindex] = 1.0/(1.0+exp((7.488-v)/5.98));

/*
		// ik1 
		ik1.Tk1ss[vindex] = 1.0/(1.0+exp(-(v+2.5538*var.ko+144.59)/(1.5692*var.ko+3.8115)));
		ik1.Ttauk1[vindex] = 122.2/(exp(-(v+127.2)/20.36)+exp((v+236.8)/69.33));
		ik1.Trk1[vindex] = 1.0/(1.0+exp((v+105.8-2.6*var.ko)/9.493));

		// inaca
		var.Thca[vindex] = exp(var.qca*v/var.RTonF);
		var.Thna[vindex] = exp(var.qna*v/var.RTonF);

		// inak 
		inak.Tknai[vindex] = inak.ko_nai*exp((inak.delta*v*F)/(3.0*R*T));
		inak.Tknao[vindex] = inak.ko_nao*exp(((1.0-inak.delta)*v*F)/(3.0*R*T));

		// ikb
		ikb.Txkb[vindex] = 1.0/(1.0+exp(-(v-14.48)/18.34));

		// icab
		icab.Texp[vindex] = exp(v/var.RTon2F);

		// inab
		inab.Texp[vindex] = exp(v/var.RTonF);

*/
	}

}
