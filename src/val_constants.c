/* produced by Tsumoto. K 2008.10.27 */

#include "syspara.h"

void val_consts(FILE *fp1)
{
	int i,w;
	double v_old,dvdt,dvdt_new;


	// Cell Geometry */
		var.length = 0.01;	// Length of the cell (cm)
		var.a = 0.0011;		// Radius of the cell (cm)
		var.vcell = 1000*M_PI*var.a*var.a*var.length; // Cell Volume:3.801e-5 (uL)
		printf("vcellm=%e\n",var.vcell);
		var.ageo = 2*M_PI*var.a*(var.a+var.length);  // geometric membrane area: 7.671e-5 (um^2)
		var.acap = var.ageo*2;          // Capacitive membrane area: 1.534e-4 cm^2 (cm^2)
		printf("Cm=%e\n",var.acap);
		var.vmyo = var.vcell*0.68;      // Myoplasm volume (uL) = 68% for Cell volume
		var.vmito = var.vcell*0.26;     // Mitochondria volume (uL) = 26% for cell volume
		var.vsr = var.vcell*0.06;       // SR volume (uL)
		var.vnsr = var.vcell*0.0552;    // NSR volume (uL)
		var.vjsr = var.vcell*0.0048;    // JSR volume (uL)
		var.vcleft = var.vcell*0.12/0.88;  // Cleft volume (uL)
		printf("const=%e\n",var.acap/var.vmyo/F);
		//exit(0);

	// Max conductanve (constant)
		// fast sodium current
		var.gna_max = 16;    // mS/uF

		// L-type calcium current
		var.pcana = 6.75E-7; // cm/s
		var.pcak = 1.93E-7; // cm/s
		var.pca = 5.4E-4; // cm/s
		var.km_ca = 0.6E-3; // mmol/L
		//var.km_ca = 0.0006;     // Half-saturation concentration of Ca channel (mM)

		// time-dependent potassium current
		var.gkx_max = 0.282;   // mS/uF
		
		// inward rectifier potassium current
		var.gk1_max = 0.75;  // mS/uF (control)
		
		// plateau potassium current
		var.gkp = 0.0183;  // mS/uF

		// Sodium-Calcium Exchanger V-S
		var.knaca = 2000.0; //(pA/pF) // control 
		var.ksat = 0.1;   // the saturation factor of inaca at very negative potentials
		var.eta = 0.35;   // position of the energy barrier controling voltage dependence of Inaca
		var.kmna = 87.5; //(mmol/L) 
		var.kmca = 1.38; //(mmol/L)

		// Sodium-Potassium Pump
		var.ibarnak = 1.5;   // Max. current through Na-K pump (uA/uF) (control)
		var.kmnai = 10;    // Half-saturation concentration of NaK pump (mmol/L)
		var.kmko = 1.5;    // Half-saturation concentration of NaK pump (mmol/L)

		// nonspecific Ca-activated current
		var.pnsca = 1.75E-7;	// (cm/s)
		var.Km_nsca = 1.2E-3;	// (mmol/L)

		// Sarcolemmal Ca pump
		var.pca_max = 1.15;	// (uA/uF)
		var.kmpca = 0.5E-3;	// (mmol/L)

		// Ca Background Current 
		var.gcab = 0.003016;  // Max. conductance of Ca background (mS/uF) (control)

	// Na Background Current 
		var.gnab = 0.00141;  // Max. conductance of Na background (mS/uF) (control)

	// NSR Ca Ion Concentration Changes 
		//var.kmup = 0.00092;   // Half-saturation concentration of iup (mM/L)
		var.kmup = 0.92E-3;   // Half-saturation concentration of iup (mmol/L)
		var.iupbar = 0.005; // Max. current through iup channel (mmol/L/ms) (control)
		var.kleak = 0.005/15.0; // Rate constant of Ca leakage from NSR (ms^-1) (control)

	// Translocation of Ca Ions from NSR to JSR
		var.tautr = 180;      // Time constant of Ca transfer from NSR to JSR (ms)
	
	// Ca buffers in the Myoplasm: Troponin (TRPN) and Calmodulin (CMDN) 
		var.cmdnbar = 0.050;   // Max. [Ca] buffered in CMDN (mM/L)
		var.trpnbar = 0.070;   // Max. [Ca] buffered in TRPN (mM/L)
		var.kmcmdn = 0.00238;  // half-saturation concentration of CMDN (mM/L)
		var.kmtrpn = 0.0005;   // half-saturation concentration of TRPN (mM/L)

	// Ca buffers in JSR: Calsequestrin (CSQN)
		var.csqnbar = 10.0;     // Max. [Ca] buffered in CSQN (mM/L)
		var.kmcsqn = 0.8;     // half-saturation concentration of CSQN (mM/L)
		
		//var.trpn = 0.0143923;
		//var.cmdn = 0.00257849;
		//var.csqn = 6.97978;
	
	// JSR Ca Ion Concentration Changes 
		var.tauon = 2.0;        // Time constant of activation of Ca release from JSR (ms)
		var.tauoff = 2.0;       // Time constant of deactivation of Ca release from JSR (ms)
		var.csqnth = 0.7;    	// orginal value. Threshold for release of Ca from CSQN due to JSR ovreload (mM)
		var.gmaxrel = 60;    	// Max. rate constant of Ca release from JSR due to overload (ms^-1)
		var.gmaxrel_jsr_ol = 0; // Rate constant of Ca release from JSR due to overload (ms^-1)
		var.swspontan = 0;    	// switch of spontaneous release
		var.km_rel = 0.8E-3;	// (mmol/L)

	// Extracellular Ion Concentrations 
		var.nao=140;
		var.ko=5.4;
		var.cao=1.8;

	// Another parameter initial setting
		var.t_cicr = 25;	//(ms)
		var.t_jsr_ol = 25;	//(ms)
		var.timer = 25;		//(ms)
		var.boolien = 0;
		var.Cainflux2 = 0.0;	//(mmol/L)
		var.Cainflux_th = 0.18E-3;	// (mmol/L)

		printf("Istim=%lf\n",var.Istim_base);
		printf("csqn=%lf\n",var.csqn);

}

