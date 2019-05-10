#include <iostream>
#include <cstring>
#include <math.h>
#include <cstdlib>
#include "Librerias.h"

using namespace std;

void MovParabolico(float VelocidadIni, float PosIniY, double Angulo, float Tiempo){
	double X=0, Y=0, VX=0, VY=0, t=0, a=0, b=0, c=0, t2=0;
	Angulo = Angulo *(M_PI/180);
	VX = VelocidadIni*cos(Angulo);
	for(int i=0;i<=Tiempo;i+=2)
	{
		if((i+2)>Tiempo){
			i=Tiempo;
		}
		if(PosIniY==0)
		{
			t = ((VelocidadIni*sin(Angulo))/9.8)*2;	
		}
		else
		{
			a = -0.5*9.8;
			b = VelocidadIni*sin(Angulo);
			c = PosIniY;
			t = (-(b)+sqrt((pow(b,2))-(4*a*c)))/(2*a);
			t2 = (-(b)-sqrt((pow(b,2))-(4*a*c)))/(2*a);
			if(t<0)
			{
				t=t2;
			}
		}
		VY = (VelocidadIni*sin(Angulo))-(9.8*i);
		X = VX*i;
		Y = PosIniY+(VelocidadIni*sin(Angulo)*i)-(0.5*9.8*pow(i,2));
		if (i>=t)
		{
			X = VX*t;
			VX=0;
			VY=0;
			Y=0;
		}
		cout<<"\nLa velocidad (eje X) en el tiempo "<<i<<" es: "<<VX;
		cout<<"\nLa posicion de la particula (eje X) en el tiempo "<<i<<" es: "<<X<<endl;
		cout<<"La velocidad (eje Y) en el tiempo "<<i<<" es: "<<VY<<endl;
		cout<<"La posicion de la particula (eje Y) en el tiempo "<<i<<" es: "<<Y<<endl;
		VX = VelocidadIni*cos(Angulo);
    }
}

void ChoqueElastico(float VelIniM1, float VelIniM2, float Masa1, float Masa2){
	float velfinM1 = 0, velfinM2 = 0;
	if((VelIniM1!=0) && (VelIniM2==0))
	{
		velfinM1 = ((Masa1-Masa2)/(Masa1+Masa2))*VelIniM1;
		velfinM2 = ((2*Masa1)/(Masa1+Masa2))*VelIniM1;
		cout<<"\nLa velocidad final de la masa 1 es: "<<velfinM1<<endl;
		cout<<"La velocidad final de la masa 2 es: "<<velfinM2<<endl;
	}
	if((VelIniM1==0) && (VelIniM2!=0))
	{
		velfinM2 = ((Masa2-Masa1)/(Masa1+Masa2))*VelIniM2;
		velfinM1 = ((2*Masa2)/(Masa1+Masa2))*VelIniM2;
		cout<<"\nLa velocidad final de la masa 1 es: "<<velfinM1<<endl;
		cout<<"La velocidad final de la masa 2 es: "<<velfinM2<<endl;
	}
	if((VelIniM1!=0) && (VelIniM2!=0))
	{
		velfinM1 = (((Masa1-Masa2)/(Masa1+Masa2))*VelIniM1) + (((2*Masa2)/(Masa1+Masa2))*VelIniM2);
		velfinM2 = (((2*Masa1)/(Masa1+Masa2))*VelIniM1) - (((Masa1-Masa2)/(Masa1+Masa2))*VelIniM2);
		cout<<"\nLa velocidad final de la masa 1 es: "<<velfinM1<<endl;
		cout<<"La velocidad final de la masa 2 es: "<<velfinM2<<endl;
	} 
}

void Rebotes(float h, float m){
	float U1=0, V1=0, x=h;
	int n=1;
	for(float i=0; i<1;i+=0.2)
	{
		cout<<"\nLa cantidad de rebotes con un coeficiente de restitucion de "<<i<<" son: "<<endl;
		if(i==0)
		{
			cout<<"Hay 0 rebotes ya que el choque es perfectamnte inelastico"<<endl;
			
		}
		else
		{
			while(h>0){
				U1=sqrt(2*9.8*h);
				V1=i*U1;
				h=(pow(i,(2*n)))*h;
				if(h>0)
				{
					cout<<"\nEl rebote "<<n<<" alcanzo una altura de: "<<h<<endl;
					n++;
				}
			}
		}
		h=x;
		n=1;
	}
}
