/* produced by Tsumoto. K 2008.10.27 */

#include <string.h>
#include <stdlib.h>
#include "syspara.h"

int mode = 1;
int P = 2;
FILE *fopen(), *fpin, *fp0, *fp1, *fp2, *fp3, *fp4, *fp5, *fp6;
int beats = 30;

typedef double Number;

main(argc,argv)
int argc;
char **argv;
{
	int i,k,w,count;
	int ii=0;
	double x[NN];
	double t = 0.0;
	double time;
	double h;
	double t_stok;
	double base_vm;
	char *tmpname;
	char cmd[BUFSIZ];
	double tend;
	double v_old,dvdt,dmax_time,pre_dvdtmax;


/* Action Potential Duration and Max. Info */
	double *vmax ; // Max. Voltage (mV)
	double *dvdtmax ; // Max. dv/dt (mV/ms)
	double *dvdtmax2 ; // Max. dv/dt (mV/ms)
	double *apd; // Action Potential Duration
	double *toneapd; // Time of dv/dt Max.
	double *ttwoapd; // Time of 90% Repolarization
	double *rmbp; // Resting Membrane Potential
	double *nair; // Intracellular Na At Rest
	double *cair; // Intracellular Ca At Rest
	double *kir ; // Intracellular K At Rest
	double caimax [beats] ; // Peak Intracellular Ca
	
	vmax=(Number *)calloc(beats,sizeof(Number));
	dvdtmax=(Number *)calloc(beats,sizeof(Number));
	dvdtmax2=(Number *)calloc(beats,sizeof(Number));
	apd=(Number *)calloc(beats,sizeof(Number));
	toneapd=(Number *)calloc(beats,sizeof(Number));
	ttwoapd=(Number *)calloc(beats,sizeof(Number));
	rmbp=(Number *)calloc(beats,sizeof(Number));
	nair=(Number *)calloc(beats,sizeof(Number));
	cair=(Number *)calloc(beats,sizeof(Number));
	kir=(Number *)calloc(beats,sizeof(Number));
	if(vmax==NULL || dvdtmax==NULL || apd==NULL || toneapd==NULL || ttwoapd==NULL || rmbp==NULL || nair==NULL || cair==NULL || kir==NULL) exit(1);

	tmpname = "temp";

	sprintf(cmd, "/usr/bin/cpp -P %s > %s", argv[1],tmpname);
	if(system(cmd) == -1){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if((fpin=fopen(tmpname,"r"))==NULL){
		fprintf(stderr,"cannot open %s\n",argv[1]);
		exit(1);
	}
	if ((fp1 = fopen("para.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp2 = fopen("data.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp3 = fopen("ndata.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp4 = fopen("current.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp5 = fopen("ikr_act.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}
	if ((fp6 = fopen("dvdt.out","w")) == NULL){
		printf("Can't open File\n");
		exit(1);
	}

	input_para(fpin);

	if (var.write){
		if ((fp0 = fopen(argv[2],"w"))==NULL){
			fprintf(stderr, "%s cannot open.\n",argv[2]);
			exit(-1);
		}
	}

	xhplot(WINDOW, 700.0, 700.0, WHITE);
	xhplot(DIRECT, 0.0, 0.0, WHITE);

	for (ii = 0; ii < var.datas; ii++){
		long j;
		time = 0.0;
		tend = var.tend[ii];
		for (i = 0; i < NN; i++){ 
			x[i] = var.x0[ii][i];
		}
		h = 1.0 / var.m;
		h *= var.tsign[ii];
		xddp.line_wid = var.line_wid[ii];
		xhplot(LINEATT,0,0,WHITE);

		// initial values input.
		val_consts(fp1);
		printf("exit consts\n");

		// initial values input.
		initial_mem();
		printf("exit memory initialization\n");

		// Tablize exp functions.       
		printf("start tablization\n");
		make_ExpTable();
		printf("finished tablization\n");

		// Initialization time
		var.dt = h;
		var.beat = 0;
		var.dvdt = 0;
		var.max_flag = 0;

		var.l = var.BCL;

		while(count < 30){
			eventloop(fp1,&mode,&P,x);
			
			var.boolien = 0;
			var.Cainflux2 = 0.0;

			for (j=0; j<(int)var.m*var.l;j++){
				t = h*(double)j;
				if ( time-(var.BCL*var.beat+10.0) >= 0.0 && time-(var.BCL*var.beat+10.0) < h ){
					apd[var.beat] =0; toneapd[var.beat] =0; ttwoapd[var.beat] =0; 
					rmbp[var.beat] =x[0]; nair[var.beat] = x[7]; kir[var.beat] = x[8]; cair[var.beat] = x[9]; caimax[var.beat] = x[9];
					vmax[var.beat] = 0.0; dvdtmax[var.beat] = -100.0; dvdtmax2[var.beat] = -100.0;
					printf("%d apd=%lf rest=%lf\n",var.beat,apd[var.beat-1],rmbp[var.beat]);
					printf("time=%lf,Istim=%lf\n",time,var.Istim);
					printf("dvdtmax[%d]=%lf\n",var.beat,dvdtmax[var.beat]);
				}
				
				if (time-(var.BCL*var.beat+10.0) >= 0.0 && time-(var.BCL*var.beat+10.0) < 0.5 ){
					var.Istim = var.Istim_base;
					if(time-(var.BCL*var.beat+10.0) <= h) printf("time=%lf,Istim=%lf\n",time,var.Istim);
				} else {
					var.Istim = 0;
				}

				v_old = x[0];
				eular(NN,h,x,t);
				var.dvdt = (x[0]-v_old)/h;
				//for(k=0;k<NN;k++){ printf("x[%d]=%lf ",k,x[k]);}printf("\n");

				if (x[0] > vmax[var.beat])
					vmax[var.beat] = x[0];
				if (x[12] > caimax[var.beat] )
					caimax[var.beat] = x[12];
				if (var.dvdt > dvdtmax[var.beat] ){
					dvdtmax[var.beat] = var.dvdt;
					toneapd[var.beat] = time;
					dmax_time = time;
				}
				if(time-(var.BCL*var.beat+10.0) >= 0.0 && time - dmax_time == 0.0 && dvdtmax[var.beat] > 100.0){
					if(fabs(pre_dvdtmax-dvdtmax[var.beat])<1E-1){
						var.boolien = 1;
						var.timer = 0.0;
					}
				}
				if (time-(var.BCL*var.beat+100.0) >= 0.0 ){
					if (var.dvdt > dvdtmax2[var.beat] ){
						dvdtmax2[var.beat] = var.dvdt;
					}
				}
				//base_vm=vmax[var.beat]-rmbp[var.beat];
				base_vm=vmax[var.beat]-rmbp[var.beat];
				if (var.dvdt < 0 && x[0] >= (vmax[var.beat] -0.9*base_vm ) ) ttwoapd[var.beat] = time;
				
				if (var.pflag) orbit(&mode,x,var.dvdt);

				//if (time>= (beats-3)*var.BCL && time < beats*var.BCL){
				if (count>25){
				//if (time>= 0.0){
					fprintf(fp2,"%lf %e %lf %lf %lf %lf %lf %lf %e %e %lf %lf\n",
							time,x[0],var.Cainflux2,var.I_all,x[7],x[8],x[9],var.t_cicr,var.Irel_cicr,var.Irel_jsr_ol,x[10],x[11]);
					fprintf(fp4,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
							time,x[0],var.ina,ical.ca,ical.na,ical.k,var.ikx,var.ik1,var.inak,var.incx,var.ipca,var.inab,var.icab);
					fprintf(fp6,"%lf %lf\n",time,var.dvdt);
					//printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
					//		time,x[0],var.ina,ical.ca,ical.na,ical.k,var.ikx,var.ik1,var.inak,var.incx,var.ipca,var.inab,var.icab);
				}
				
				time += h;
				var.clock += h;
				var.t_cicr += h;
				var.t_jsr_ol += h;
				var.timer += h;
				pre_dvdtmax=dvdtmax[var.beat];

			} // end j-loop

			draw_p(&mode,P,x,var.dvdt);
			mouse(&mode,x,var.dvdt);
			if (fabs(time) > tend &&  tend != 0.0) break;
			var.beat++;
			count++;

		} // end while loop

		fclose(fp1);
		fclose(fp2);
		fclose(fp3);
		fclose(fp4);
		fclose(fp5);
		fclose(fp6);
		free(vmax);free(dvdtmax);free(dvdtmax2);free(apd);free(toneapd);free(ttwoapd);
		free(rmbp);free(nair);free(cair);free(kir);
		closed_mem();

	} // end ii-loop
} // end main

