#include "syspara.h"

void comp_rev_potential(double x[])
{

	var.Ena=((R*T)/F)*log(var.nao/x[7]); // [Na]i = x[7],[Na]o = var.nao
	var.Ekx = (R*T/F)*log((var.ko+0.01833*var.nao)/(x[8]+0.01833*x[7]));
	var.Ek = (R*T/F)*log(var.ko/x[8]);
	var.Eca=(R*T)/(2*F)*log(var.cao/x[9]);
	var.Ensca = ((R*T)/F)*log((var.ko + var.nao)/(x[7] + x[8]));

}

void comp_ina(double x[])
{

	int iV=0;
	double V1,V2,d1,d2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ina.mss = ina.Tmss[iV]*d2 + ina.Tmss[iV+1]*d1;
	ina.taum = ina.Ttaum[iV]*d2 + ina.Ttaum[iV+1]*d1;
	ina.hss = ina.Thss[iV]*d2 + ina.Thss[iV+1]*d1;
	ina.tauh = ina.Ttauh[iV]*d2 + ina.Ttauh[iV+1]*d1;
	ina.jss = ina.Tjss[iV]*d2 + ina.Tjss[iV+1]*d1;
	ina.tauj = ina.Ttauj[iV]*d2 + ina.Ttauj[iV+1]*d1;

	var.ina = var.gna_max*x[1]*x[1]*x[1]*x[2]*x[3]*(x[0]-var.Ena);
	//printf("Vm=%lf ina=%lf\n",x[0],var.ina);
}


void comp_ical(double x[])
{

	int iV=0;
	double V1,V2,d1,d2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ical.dss = ical.Tdss[iV]*d2 + ical.Tdss[iV+1]*d1;
	ical.taud = ical.Ttaud[iV]*d2 + ical.Ttaud[iV+1]*d1;
	ical.fss = ical.Tfss[iV]*d2 + ical.Tfss[iV+1]*d1;
	ical.tauf = ical.Ttauf[iV]*d2 + ical.Ttauf[iV+1]*d1;
	
	var.fca = 1.0/(1.0+(x[9]/var.km_ca)*(x[9]/var.km_ca));
	
	ical.namax=var.pcana*1.0*(x[0]*F*F/(R*T))*(0.75*x[7]*exp(x[0]*F/(R*T)) - 0.75*var.nao)/(exp(x[0]*F/(R*T)) - 1.0);
	ical.kmax=var.pcak*1.0*(x[0]*F*F/(R*T))*(0.75*x[8]*exp(x[0]*F/(R*T)) - 0.75*var.ko)/(exp(x[0]*F/(R*T)) - 1.0);
	ical.camax=var.pca*4.0*(x[0]*F*F/(R*T))*(1.0*x[9]*exp(2.0*x[0]*F/(R*T)) - 0.341*var.cao)/(exp(2.0*x[0]*F/(R*T))-1.0);

	ical.na=ical.namax*x[4]*x[5]*var.fca;
	ical.k=ical.kmax*x[4]*x[5]*var.fca;
	ical.ca=ical.camax*x[4]*x[5]*var.fca;

	ical.total = ical.na + ical.k + ical.ca;

	//printf("Vm=%lf ical=%lf\n",x[0],ical.ca);
}

// Time-dependent Potassium Current 
void comp_ikx (double x[])
{
	
	int iV=0;
	double V1,V2,d1,d2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikx.xss = ikx.Txss[iV]*d2 + ikx.Txss[iV+1]*d1;
	ikx.taux = ikx.Ttaux[iV]*d2 + ikx.Ttaux[iV+1]*d1;
	ikx.rik = ikx.Trik[iV]*d2 + ikx.Trik[iV+1]*d1;
	
	var.ikx =  var.gkx_max*sqrt(var.ko/5.4)*ikx.rik*x[6]*x[6]*(x[0]-var.Ekx);

}

void comp_ik1 (double x[])
{

	var.Ek1 = var.Ek;

	var.ak1 = 1.02/(1.0+exp(0.2385*(x[0]-var.Ek1-59.215)));
	var.bk1 = (0.49124*exp(0.08032*(x[0]-var.Ek1+5.476))+exp(0.06175*(x[0]-var.Ek1-594.31)))/(1.0+exp(-0.5143*(x[0]-var.Ek1+4.753)));
	var.k1_inf = var.ak1/(var.ak1+var.bk1);

	var.ik1 = var.gk1_max*sqrt(var.ko/5.4)*var.k1_inf*(x[0]-var.Ek1);

}

void comp_ikp (double x[])
{

	int iV=0;
	double V1,V2,d1,d2;

	V1 = (x[0]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikp.xkp = ikp.Txkp[iV]*d2 + ikp.Txkp[iV+1]*d1;
	
	var.ikp = var.gkp*ikp.xkp*(x[0]-var.Ek);

}

// Sodium-Calcium Exchanger V-S
void comp_incx (double x[])
{

	double para;
	double u1,u2,u3,u4;

	para = F/(R*T);

	u1 = 1.0/(var.kmna*var.kmna*var.kmna + var.nao*var.nao*var.nao);
	u2 = 1.0/(var.kmca + var.cao);
	u3 = 1.0/(1.0+var.ksat*exp((var.eta-1.0)*x[0]*para));
	u4 = exp(var.eta*x[0]*para)*x[7]*x[7]*x[7]*var.cao - exp((var.eta-1.0)*x[0]*para)*var.nao*var.nao*var.nao*x[9];

    var.incx = var.knaca*u1*u2*u3*u4;

}

// Sodium-Potassium Pump
void comp_inak (double x[])
{

	double para;

	para = F/(R*T);

	var.sigma = (exp(var.nao/67.3)-1.0)/7.0;
	var.fnak = 1.0/(1.0+0.1245*exp(-0.1*x[0]*para)+0.0365*var.sigma*exp(-x[0]*para));

	var.inak = var.ibarnak*var.fnak*1.0/(1.0+pow(var.kmnai/x[7],1.5))*var.ko/(var.ko+var.kmko);

}

// Nonspecific Ca-activated current
void comp_insca (double x[])
{

	double para;

	para = F/(R*T);
	var.vns = x[0] - var.Ensca;

	var.pnsna = var.pnsca*(var.vns*F*para)*(0.75*x[7]*exp(var.vns*para)-0.75*var.nao)/(exp(var.vns*para) - 1.0);
    var.pnsk  = var.pnsca*(var.vns*F*para)*(0.75*x[8]*exp(var.vns*para)-0.75*var.ko)/(exp(var.vns*para) - 1.0);

	var.insna = var.pnsna/(1.0+(var.Km_nsca/x[9])*(var.Km_nsca/x[9])*(var.Km_nsca/x[9]));
	var.insk  = var.pnsk/(1.0+(var.Km_nsca/x[9])*(var.Km_nsca/x[9])*(var.Km_nsca/x[9]));
	var.insca = var.insna + var.insk;

}

// Sarcolemmal Ca pump
//double pca_max = 1.15; // (uA/uF)
//double kmpca=0.5E-3;	// (mM/L)
void comp_ipca (double x[])
{

	var.ipca = var.pca_max*x[9]/(var.kmpca + x[9]);

}

// Ca Background Current 
// gcab = 0.003016;  // Max. conductance of Ca background (mS/uF)
// Eca;  // Nernst potential for Ca (mV)

void comp_icab (double x[])
{

	var.icab = var.gcab*(x[0]-var.Eca);

}

// Na Background Current 
// gnab = 0.00141;  // Max. conductance of Ca background (mS/uF)
// Ena;  // Nernst potential for Ca (mV)

void comp_inab (double x[])
{

	var.inab = var.gnab*(x[0]-var.Ena);

}

void conc_nsr (double x[])

// NSR Ca Ion Concentration Changes 
// iup;      // Ca uptake from myo. to NSR (mM/L/ms)
// ileak;    // Ca leakage from NSR to myo. (mM/L/ms)
// itr;      // Translocation current of Ca ions from NSR to JSR (mM/L/ms)
// x[9] = [Ca]i (cai)
// x[10] = [Ca]_NSR (nsr)
// x[11] = [Ca]_JSR (jsr)

{

	var.iup = var.iupbar*x[9]/(x[9]+var.kmup);
	var.ileak = var.kleak*x[10];
	var.itr = (x[10]-x[11])/var.tautr; 
}

void conc_cai (double x[])

// Myoplasmic Ca Ion Concentration Changes 
// catotal; // Total myoplasmic Ca concentration (mM)
// cmdnbar = 0.050;   // Max. [Ca] buffered in CMDN (mM/L)
// trpnbar = 0.070;   // Max. [Ca] buffered in TRPN (mM/L)
// kmcmdn = 0.00238;  // Equalibrium constant of buffering for CMDN (mM/L)
// kmtrpn = 0.0005;   // Equalibrium constant of buffering for TRPN (mM/L)

{

	double u1,u2,u3,u4,u5,u6;

	u1 = (var.kmcmdn+x[9])*(var.kmcmdn+x[9]);
	u2 = (var.kmtrpn+x[9])*(var.kmtrpn+x[9]);
	u3 = var.cmdnbar*var.kmcmdn/u1;
	u4 = var.trpnbar*var.kmtrpn/u2;

	u5 = (var.kmcsqn+x[11])*(var.kmcsqn+x[11]);
	u6 = var.csqnbar*var.kmcsqn/u5;

	var.cai_bufc = 1.0/(1.0+u3+u4);

	var.ca_jsr_bufc = 1.0/(1.0+u6);

	var.csqn = var.csqnbar*(x[11]/(x[11]+var.kmcsqn));

	//if(var.csqn>8.0){
	//	printf("ca_jsr=%e,kmcaqn=%lf,buf_csqn=%lf\n",x[14],var.kmcsqn,var.csqn);
	//}
	
}

void conc_jsr (double x[])

// JSR Ca Ion Concentration Changes 
// tauon = 2;        	// Time constant of activation of Ca release from JSR (ms)
// tauoff = 2;       	// Time constant of deactivation of Ca release from JSR (ms)
// magrel;           	// Magnitude of Ca release
// irelcicr;         	// Ca release from JSR to myo. due to CICR (mM/ms)
// csqnth = 0.7;    	// Threshold for release of Ca from CSQN due to JSR ovreload (mM)
// csqnth = 7.0;    	// Threshold for release of Ca from CSQN due to JSR ovreload (mM)
// gmaxrel = 22;    	// Max. rate constant of Ca release from JSR due to overload (ms^-1)
// grelbarjsrol;     	// Rate constant of Ca release from JSR due to overload (ms^-1)
// greljsrol;        	// Rate constant of Ca release from JSR due to CICR (ms^-1)
// ireljsrol;        	// Ca release from JSR to myo. due to JSR overload (mM/ms)
// swspontan = 0;    	// switch of spontaneous release
// csqnbar = 10;     	// Max. [Ca] buffered in CSQN (mM)
// kmcsqn = 0.8;     	// Equalibrium constant of buffering for CSQN (mM)
// on;               	// Time constant of activation of Ca release from JSR (ms)
// off;              	// Time constant of deactivation of Ca release from JSR (ms)
// dICa;          		// Rate of change of Ca entry
// dICa_new;       		// New rate of change of Ca entry
// ICa_total_old;        // Old rate of change of Ca entry
// x[9] = [Ca]i (cai)
// x[10] = [Ca]_NSR (nsr)
// x[11] = [Ca]_JSR (jsr)
// tcicr = t=0 at time of CICR (ms)
// tjsrol = t=0 at time of JSR overload (ms)

{

	double magrel;
	double grel,greljsrol;
	double t_stok,para;

	para = var.acap/(2.0*var.vmyo*F);
	
	if( var.boolien == 1 ){
		var.timer = 0.0;
		var.boolien = 2;
		var.Cainflux2 = 0.0;
	}

	if(var.boolien == 2){
		var.Cainflux2 += var.dt*(-(ical.ca + var.icab + var.ipca - 2.0*var.incx)*para);
		//printf("%d Cainf2=%e timer=%lf time=%lf vm=%lf dt=%lf\n",var.boolien,var.Cainflux2,var.timer,var.clock,x[0],var.dt);
	}
	
	if(fabs(var.timer-2.0)<0.00001 && var.boolien==2){
		if(var.Cainflux2 > var.Cainflux_th){
			var.boolien = 3;
			var.t_cicr = 0.0;
			var.magrel = (var.Cainflux2-var.Cainflux_th)/(var.km_rel+var.Cainflux2-var.Cainflux_th);
			printf("start CICR. magrel=%e, tcicr=%lf, Grelmax=%lf\n",magrel,var.t_cicr,var.gmaxrel);
		} else {
			var.magrel += 0.0;
		}
		//printf("%d Cainf2=%e timer=%lf time=%lf vm=%lf dt=%lf\n",var.boolien,var.Cainflux2,var.timer,var.clock,x[0],var.dt);
	}
	
	grel = var.gmaxrel*var.magrel*(1.0-exp(-var.t_cicr/var.tauon))*exp(-var.t_cicr/var.tauoff);
	var.Irel_cicr = grel*(x[11]-x[9]);

	var.csqn = var.csqnbar*(x[11]/(x[11]+var.kmcsqn));

	if(var.csqn >= var.csqnth ){
		var.t_jsr_ol = 0.0;
		var.gmaxrel_jsr_ol = 4.0;
	} else {
		var.gmaxrel_jsr_ol = 0.0;
	}
	
	greljsrol = var.gmaxrel_jsr_ol*(1.0-exp(-var.t_jsr_ol/var.tauon))*exp(-var.t_jsr_ol/var.tauoff);
	
	var.Irel_jsr_ol = greljsrol*(x[14]-x[12]);
	

}

