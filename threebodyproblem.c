#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define w 2*M_PI

typedef struct {
double x;
double y;
double vx;
double vy;
} nuovotipo;

typedef struct {
double U;
double K;
double E;
double J;
double v2;
} nuovotipo2;

void rungekutta2(FILE *fp,nuovotipo old,double Mb,double dt,double xa,double xb,double tmax);
void errore(FILE* fq,nuovotipo old,double Mb,double dt,double xa,double xb,double tmax,nuovotipo ci,nuovotipo2 di);
double ax(double Mb,nuovotipo z,double xa, double xb);
double ay(double Mb,nuovotipo z,double xa, double xb);
double dist(double x2,double x1,double y);
nuovotipo2 energie(double Mb, nuovotipo a,double xa,double xb);

int main(){
	/*if(argc!=7){
	fprintf(stderr, "\n ordine di inserimento:\n <Mb> <x(0)> <y(0)> <vx(0)> <vy(0)> <dt>\n",argv[0]);
	exit(1);
	}*/
	double Ma,Mb,dt,t=0;
	double xa,xb,ya,yb;
	nuovotipo c;
	nuovotipo2 d;
	nuovotipo ci;
	nuovotipo2 di;
	int i;
	
	/*Mb=strtod(argv[1],NULL);
	c.x=strtod(argv[2],NULL);
        c.y=strtod(argv[3],NULL);
	c.vx=strtod(argv[4],NULL);
	c.vy=strtod(argv[5],NULL);
	dt=strtod(argv[6],NULL);*/
	
	printf("inserisci la massa Mb \n");
        scanf("%lf",&Mb);
        printf("inserisci x(0) \n");
        scanf("%lf",&c.x);
        printf("inserisci y(0) \n");
        scanf("%lf",&c.y);
	printf("inserisci vx(0) \n");
        scanf("%lf",&c.vx);
        printf("inserisci vy(0) \n");
        scanf("%lf",&c.vy);
        printf("inserisci il passo d'integrazione\n");
  	scanf("%lf",&dt);
  	
  	Ma=1-Mb;
	xa=-Mb;
	ya=0;
	xb=Ma;
	yb=0;
  	
  	ci.x=c.x;
  	ci.y=c.y;
  	ci.vx=c.vx;
  	ci.vy=c.vy;
  	
  	di.v2=ci.vx*ci.vx+ci.vy*ci.vy;
  	di.U=(-4*M_PI*M_PI*((Ma/dist(xa,ci.x,ci.y))+(Mb/dist(xb,ci.x,ci.y))));
	di.K=0.5*di.v2;
	di.E=di.U+di.K;
	di.J=ci.vx*ci.vx+ci.vy*ci.vy+2*di.U-w*w*(ci.x*ci.x+ci.y*ci.y);
	
	FILE *fp;
 		
	FILE *fq;
	
	rungekutta2(fp,c,Mb,dt,xa,xb,10);
	errore(fq,c,Mb,0.001,xa,xb,2,ci,di);

}

void rungekutta2(FILE *fp,nuovotipo old,double Mb,double dt,double xa,double xb,double tmax){
	fp=fopen("datif.dat","w+");
	double xplus,yplus,vxplus,vyplus;
	double Ma=1-Mb;
	double t=0;
	nuovotipo new;
	nuovotipo temp;
	fprintf(fp,"#t\t      x\t         y\t    vx\t       vy\t      U\t     K\t       E\t        J\n");
	while(t<tmax){
	nuovotipo2 d=energie(Mb,old,xa,xb);
			
	xplus=old.vx*dt;
	yplus=old.vy*dt;
	
	vxplus=ax(Mb,old,xa,xb)*dt;
	vyplus=ay(Mb,old,xa,xb)*dt;
	
	temp.x=old.x+0.5*xplus;
	temp.y=old.y+0.5*yplus;
	temp.vx=old.vx+0.5*vxplus;
	temp.vy=old.vy+0.5*vyplus;
	
	new.x=old.x+(old.vx+0.5*vxplus)*dt;
	new.y=old.y+(old.vy+0.5*vyplus)*dt;
	
	new.vx=old.vx+ax(Mb,temp,xa,xb)*dt;
	new.vy=old.vy+ay(Mb,temp,xa,xb)*dt;
	
fprintf(fp,"%lf   %lf   %lf    %lf    %lf    %lf    %lf    %lf    %lf\n",t,old.x,old.y,old.vx,old.vy,d.U,d.K,d.E,d.J);
	t=t+dt;
	old=new;
	}
	fclose(fp);
}

double ax(double Mb,nuovotipo z,double xa, double xb){
	double Ma=1-Mb;	
	double rAP3=dist(xa,z.x,z.y)*dist(xa,z.x,z.y)*dist(xa,z.x,z.y);
	double rBP3=dist(xb,z.x,z.y)*dist(xb,z.x,z.y)*dist(xb,z.x,z.y);

	return -4*M_PI*M_PI*(Ma*(z.x-xa)/rAP3+Mb*(z.x-xb)/rBP3)+2*w*z.vy+w*w*z.x;
}

double ay(double Mb,nuovotipo z,double xa, double xb){
	double Ma=1-Mb;
	double rAP3=dist(xa,z.x,z.y)*dist(xa,z.x,z.y)*dist(xa,z.x,z.y);
	double rBP3=dist(xb,z.x,z.y)*dist(xb,z.x,z.y)*dist(xb,z.x,z.y);

	return -4*M_PI*M_PI*z.y*(Ma/rAP3+Mb/rBP3)-2*w*z.vx+w*w*z.y;
}

double dist(double x2,double x1,double y){
	return sqrt((x2-x1)*(x2-x1)+y*y);
}

nuovotipo2 energie(double Mb, nuovotipo a,double xa,double xb){
	double Ma=1-Mb;
	double v2=a.vx*a.vx+a.vy*a.vy;
	nuovotipo2 new;
	
	new.U=(-4*M_PI*M_PI*((Ma/dist(xa,a.x,a.y))+(Mb/dist(xb,a.x,a.y))));
	new.K=0.5*v2;
	new.E=new.K+new.U;
	new.J=a.vx*a.vx+a.vy*a.vy+2*new.U-w*w*(a.x*a.x+a.y*a.y);
	
	return new;
}


void errore(FILE* fq,nuovotipo old,double Mb,double dt,double xa,double xb,double tmax,nuovotipo ci,nuovotipo2 di){
	fq=fopen("erroref.dat","w+");
	fprintf(fq,"#t\t      dt\t         J\t    J0\t      |(J-J0)/J0|\n");
	double xplus,yplus,vxplus,vyplus;
	double Ma=1-Mb;
	double t=0;
	nuovotipo new;
	nuovotipo temp;

	
do{	
	old=ci;
	t=0;
	while(t<tmax){
	nuovotipo2 d=energie(Mb,old,xa,xb);
			
	xplus=old.vx*dt;
	yplus=old.vy*dt;
	
	vxplus=ax(Mb,old,xa,xb)*dt;
	vyplus=ay(Mb,old,xa,xb)*dt;
	
	temp.x=old.x+0.5*xplus;
	temp.y=old.y+0.5*yplus;
	temp.vx=old.vx+0.5*vxplus;
	temp.vy=old.vy+0.5*vyplus;
	
	new.x=old.x+(old.vx+0.5*vxplus)*dt;
	new.y=old.y+(old.vy+0.5*vyplus)*dt;
	
	new.vx=old.vx+ax(Mb,temp,xa,xb)*dt;
	new.vy=old.vy+ay(Mb,temp,xa,xb)*dt;
	
	t=t+dt;
	old=new;
	if(t>2-0.5*dt && t<2+0.5*dt){
    	fprintf(fq,"  %lf          %.14lf        %lf    %lf   %.20lf\n",t,dt,d.J,di.J,fabs(((d.J-di.J)/di.J)));
		}
	}
	dt=dt/2;
	}while(dt>0.0000009);

fclose(fq);
}

