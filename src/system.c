#include "syspara.h"

void function(double x[],double f[],double t)
{
   	double f0,k1,k2,volrat_nsr,volrat_jsr; 
	
	k1= var.acap/var.vmyo/F;
	k2= var.acap/var.vmyo/F/2.0;
	volrat_nsr = var.vnsr/var.vmyo;
	volrat_jsr = var.vjsr/var.vmyo;

	comp_rev_potential(x);
	//printf("rev\n");
	comp_ina(x);
	//printf("ina\n");
	comp_ical(x);
	//printf("ical\n");
	comp_ikx(x);
	//printf("ikx\n");
	comp_ik1(x);
	//printf("ik1\n");
	comp_ikp(x);
	//printf("ikp\n");
	comp_incx(x);
	//printf("incx\n");
	comp_inak(x);
	//printf("inak\n");
	comp_insca(x);
	//printf("insca\n");
	comp_ipca(x);
	//printf("ipca\n");
	comp_icab(x);
	//printf("icab\n");
	comp_inab(x);
	//printf("inab\n");
	conc_nsr(x);
	//printf("nsr\n");
	conc_cai(x);
	//printf("cai\n");
	conc_jsr(x);
	//printf("jsr\n");

	var.INa_total = var.ina + ical.na + var.inab + var.insna + 3.0*var.inak + 3.0*var.incx;
	var.IK_total = var.ikx + var.ik1 + var.ikp + ical.k + var.insk - 2.0*var.inak + var.Istim;
	var.ICa_total = ical.ca + var.ipca + var.icab - 2.0*var.incx;
	var.I_all = var.INa_total+var.IK_total+var.ICa_total;

	f[0] = -(var.INa_total+var.IK_total+var.ICa_total);
	
	f[1] = (ina.mss - x[1])/ina.taum; // m: activation for Ina
	f[2] = (ina.hss - x[2])/ina.tauh; // h: inactivation for Ina
	f[3] = (ina.jss - x[3])/ina.tauj; // j: slowly inactivation for Ina

	f[4] = (ical.dss - x[4])/ical.taud; // d: activation for L-type Ca channel
	f[5] = (ical.fss - x[5])/ical.tauf; // f: inactivation for L-type Ca channel
	
	f[6] = (ikx.xss - x[6])/ikx.taux; // x: activation for Ik 

	f[7] = -var.INa_total*k1;
	f[8] = -var.IK_total*k1;
	f[9] = var.cai_bufc*(-var.ICa_total*k2+(var.ileak-var.iup)*volrat_nsr+(var.Irel_cicr+var.Irel_jsr_ol)*volrat_jsr);
	f[10] = var.iup - var.ileak - var.itr*var.vjsr/var.vnsr;
	f[11] = var.ca_jsr_bufc*(var.itr - var.Irel_cicr - var.Irel_jsr_ol);

}


