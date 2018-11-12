
#include "syspara.h"

typedef double Number;

void initial_mem()
{

	// ina_fast
	ina.Tmss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttaum=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Thss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauh=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Tjss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauj=(Number *)calloc(VNMAX,sizeof(Number));
	if( ina.Tmss==NULL || ina.Ttaum==NULL || ina.Thss==NULL || ina.Ttauh==NULL || ina.Tjss==NULL || ina.Ttauj==NULL ) exit(1);
	
	// ical
	ical.Tdss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttaud=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Tfss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttauf=(Number *)calloc(VNMAX,sizeof(Number));
	if( ical.Tdss==NULL || ical.Ttaud==NULL || ical.Tfss==NULL || ical.Ttauf==NULL ) exit(1);

	// ikx
	ikx.Txss=(Number *)calloc(VNMAX,sizeof(Number));
	ikx.Ttaux=(Number *)calloc(VNMAX,sizeof(Number));
	ikx.Trik=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikx.Txss==NULL || ikx.Ttaux==NULL || ikx.Trik==NULL ) exit(1);
	
	// ikp
	ikp.Txkp=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikp.Txkp==NULL ) exit(1);

/*
	// ik1
	ik1.Tk1ss=(Number *)calloc(VNMAX,sizeof(Number));
	ik1.Ttauk1=(Number *)calloc(VNMAX,sizeof(Number));
	ik1.Trk1=(Number *)calloc(VNMAX,sizeof(Number));
	if( ik1.Tk1ss == NULL || ik1.Ttauk1==NULL || ik1.Trk1==NULL ) exit(1);
	// inaca
	var.Thca=(Number *)calloc(VNMAX,sizeof(Number));
	var.Thna=(Number *)calloc(VNMAX,sizeof(Number));
	if( var.Thca==NULL || var.Thna==NULL ) exit(1);
	// inak
	inak.Tknai=(Number *)calloc(VNMAX,sizeof(Number));
	inak.Tknao=(Number *)calloc(VNMAX,sizeof(Number));
	if( inak.Tknai==NULL || inak.Tknao==NULL ) exit(1);
	// ikb	
	ikb.Txkb=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikb.Txkb==NULL ) exit(1);
	// icab
	icab.Texp=(Number *)calloc(VNMAX,sizeof(Number));
	if( icab.Texp==NULL ) exit(1);
	// inab
	inab.Texp=(Number *)calloc(VNMAX,sizeof(Number));
	if( inab.Texp==NULL ) exit(1);
*/

}


void closed_mem()
{

	free(ina.Tmss); free(ina.Ttaum); 
	free(ina.Thss); free(ina.Ttauh);
	free(ina.Tjss); free(ina.Ttauj);
	free(ical.Tdss); free(ical.Ttaud); free(ical.Tfss); free(ical.Ttauf);
	free(ikx.Txss); free(ikx.Ttaux); free(ikx.Trik); 
	free(ikp.Txkp); 
	/*
	free(ik1.Tk1ss); free(ik1.Ttauk1); free(ik1.Trk1);
	free(var.Thca); free(var.Thna);
	free(inak.Tknai); free(inak.Tknao);
	free(ikb.Txkb);
	free(icab.Texp);
	free(inab.Texp);
	*/

}

