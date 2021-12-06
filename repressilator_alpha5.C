#include <iostream>
#include <random>
#include <fstream>
#include "armadillo"
using namespace arma;
using namespace std;

#include <cstdio>
#include <cmath>
#include <cstdlib>

//#define N 4096
//#define N 8192
#define N 16384
#define N1 8192

double t=0.0;
//double dt=0.001;
double dt=0.01;
double dt_init=0.1;
vec m11array(N);
vec m21array(N);
vec m31array(N);
vec m11acf(N1);
vec m21acf(N1);
vec m31acf(N1);

double m11;
double m12;
double m13;

double m21;
double m22;
double m23;

double m31;
double m32;
double m33;

double p11;
double p12;
double p13;

double p21;
double p22;
double p23;

double p31;
double p32;
double p33;

double n=2.4;


double alpha1=5.0;
double alpha2=6.0;
double alpha3=7.0;


double alpha01=0.1;
double alpha02=0.1;
double alpha03=0.1;

double beta1=1.0;
double beta2=1.0;
double beta3=1.0;
double phi=1.0;

double summ11;
double summ12;
double summ13;
double summ21;
double summ22;
double summ23;
double summ31;
double summ32;
double summ33;
double sump11;
double sump12;
double sump13;
double sump21;
double sump22;
double sump23;
double sump31;
double sump32;
double sump33;

double keisu=10.0;
double km;
double kp;
int sample=100;
vec m11shuki(sample);
vec m21shuki(sample);
vec m31shuki(sample);

int kmcounter;//km=1~1000
double sum1=0.0;
double sum2=0.0;
double sum3=0.0;
double ave1=0.0;
double ave2=0.0;
double ave3=0.0;

double alpha=0.0;


int main(int argc, char* argv[])
{
  alpha2=5.0;
  if(argc > 1){
    alpha = atof(argv[1]);
    alpha1=alpha2-alpha;
    alpha3=alpha2+alpha;
  }
  std::random_device seed_gen;
  std::default_random_engine engine(seed_gen());
  
  // 0.0以上1.0未満の値を等確率で発生させる
  std::uniform_real_distribution<> dist(0.0, 1.0);
  for(kmcounter=1; kmcounter<1000; kmcounter++)
    {
      
      km=kmcounter*0.01;
      kp=kmcounter*0.01;
      
      /*
      km=0.5;
      kp=0.5;
      */
      /*
	km=2.0;
	kp=2.0;
      */
      /*
	km=5.0;
	kp=5.0;
      */
      for(int samplecounter=0; samplecounter<sample; samplecounter++)
	{
	  t=0.0;
	  
	  m11=dist(engine);
	  m12=dist(engine);
	  m13=dist(engine);
	  m21=dist(engine);
	  m22=dist(engine);
	  m23=dist(engine);
	  m31=dist(engine);
	  m32=dist(engine);
	  m33=dist(engine);
	  p11=dist(engine);
	  p12=dist(engine);
	  p13=dist(engine);
	  p21=dist(engine);
	  p22=dist(engine);
	  p23=dist(engine);
	  p31=dist(engine);
	  p32=dist(engine);
	  p33=dist(engine);
	  
	  summ11=m11;
	  summ12=m12;
	  summ13=m13;
	  summ21=m21;
	  summ22=m22;
	  summ23=m23;
	  summ31=m31;
	  summ32=m32;
	  summ33=m33;
	  sump11=p11;
	  sump12=p12;
	  sump13=p13;
	  sump21=p21;
	  sump22=p22;
	  sump23=p23;
	  sump31=p31;
	  sump32=p32;
	  sump33=p33;
	  
	  m11shuki(samplecounter)=0.0;
	  m21shuki(samplecounter)=0.0;
	  m31shuki(samplecounter)=0.0;
	  
	  
	  //int k=0;
	  for(int i=1; i<=10000; i++)
	    {
	      t=t+dt_init;
	      //repressilator1
	      double deltam11=alpha1/(1.0+pow(p13,n))+alpha01-keisu*m11/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam12=alpha1/(1.0+pow(p11,n))+alpha01-keisu*m12/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam13=alpha1/(1.0+pow(p12,n))+alpha01-keisu*m13/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltap11=beta1*m11-keisu*p11/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap12=beta1*m12-keisu*p12/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap13=beta1*m13-keisu*p13/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      
	      
	      //repressilator2
	      double deltam21=alpha2/(1.0+pow(p23,n))+alpha02-keisu*m21/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam22=alpha2/(1.0+pow(p21,n))+alpha02-keisu*m22/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam23=alpha2/(1.0+pow(p22,n))+alpha02-keisu*m23/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltap21=beta2*m21-keisu*p21/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap22=beta2*m22-keisu*p22/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap23=beta2*m23-keisu*p23/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      
	      //repressilator3
	      double deltam31=alpha3/(1.0+pow(p33,n))+alpha03-keisu*m31/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam32=alpha3/(1.0+pow(p31,n))+alpha03-keisu*m32/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam33=alpha3/(1.0+pow(p32,n))+alpha03-keisu*m33/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltap31=beta3*m31-keisu*p31/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap32=beta3*m32-keisu*p32/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap33=beta3*m33-keisu*p33/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      
	      
	      m11=m11+deltam11*dt;
	      m12=m12+deltam12*dt;
	      m13=m13+deltam13*dt;
	      m21=m21+deltam21*dt;
	      m22=m22+deltam22*dt;
	      m23=m23+deltam23*dt;
	      m31=m31+deltam31*dt;
	      m32=m32+deltam32*dt;
	      m33=m33+deltam33*dt;
	      p11=p11+deltap11*dt;
	      p12=p12+deltap12*dt;
	      p13=p13+deltap13*dt;
	      p21=p21+deltap21*dt;
	      p22=p22+deltap22*dt;
	      p23=p23+deltap23*dt;
	      p31=p31+deltap31*dt;
	      p32=p32+deltap32*dt;
	      p33=p33+deltap33*dt;
	      
	      //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",t,m11,m12,m13,m21,m22,m23,m31,m32,m33,p11,p12,p13,p21,p22,p23,p31,p32,p33);
	    }
	  sum1=0.0;
	  sum2=0.0;
	  sum3=0.0;
	  for(int i=1; i<=N; i++)
	    {
	      t=t+dt;
	      //repressilator1
	      double deltam11=alpha1/(1.0+pow(p13,n))+alpha01-keisu*m11/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam12=alpha1/(1.0+pow(p11,n))+alpha01-keisu*m12/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam13=alpha1/(1.0+pow(p12,n))+alpha01-keisu*m13/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltap11=beta1*m11-keisu*p11/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap12=beta1*m12-keisu*p12/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap13=beta1*m13-keisu*p13/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      
	      
	      //repressilator2
	      double deltam21=alpha2/(1.0+pow(p23,n))+alpha02-keisu*m21/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam22=alpha2/(1.0+pow(p21,n))+alpha02-keisu*m22/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam23=alpha2/(1.0+pow(p22,n))+alpha02-keisu*m23/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltap21=beta2*m21-keisu*p21/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap22=beta2*m22-keisu*p22/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap23=beta2*m23-keisu*p23/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      
	      //repressilator3
	      double deltam31=alpha3/(1.0+pow(p33,n))+alpha03-keisu*m31/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam32=alpha3/(1.0+pow(p31,n))+alpha03-keisu*m32/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltam33=alpha3/(1.0+pow(p32,n))+alpha03-keisu*m33/(km+m11+m12+m13+m21+m22+m23+m31+m32+m33);
	      double deltap31=beta3*m31-keisu*p31/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap32=beta3*m32-keisu*p32/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      double deltap33=beta3*m33-keisu*p33/(kp+p11+p12+p13+p21+p22+p23+p31+p32+p33);
	      
	      
	      m11=m11+deltam11*dt;
	      m12=m12+deltam12*dt;
	      m13=m13+deltam13*dt;
	      m21=m21+deltam21*dt;
	      m22=m22+deltam22*dt;
	      m23=m23+deltam23*dt;
	      m31=m31+deltam31*dt;
	      m32=m32+deltam32*dt;
	      m33=m33+deltam33*dt;
	      p11=p11+deltap11*dt;
	      p12=p12+deltap12*dt;
	      p13=p13+deltap13*dt;
	      p21=p21+deltap21*dt;
	      p22=p22+deltap22*dt;
	      p23=p23+deltap23*dt;
	      p31=p31+deltap31*dt;
	      p32=p32+deltap32*dt;
	      p33=p33+deltap33*dt;
	      m11array(i-1)=m11;
	      m21array(i-1)=m21;
	      m31array(i-1)=m31;
	      sum1=sum1+m11;
	      sum2=sum2+m21;
	      sum3=sum3+m31;
	      
	      //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",t,m11,m12,m13,m21,m22,m23,m31,m32,m33,p11,p12,p13,p21,p22,p23,p31,p32,p33);
	      //printf("%lf %lf %lf %lf\n",t,m11,m21,m31);
	    }
	  
	  ave1=sum1/(double)N;
	  ave2=sum2/(double)N;
	  ave3=sum3/(double)N;
	  //printf("%lf %lf %lf\n",ave1,ave2,ave3);
	  for(int i=0; i<N; i++)
	    {
	      m11array(i)=m11array(i)-ave1;
	      m21array(i)=m21array(i)-ave2;
	      m31array(i)=m31array(i)-ave3;
	      //printf("%d %lf %lf %lf\n",i,m11array(i-1),m21array(i-1),m31array(i-1));
	    }
	  double phi01=0.0;
	  double phi02=0.0;
	  double phi03=0.0;
	  for(int i=0; i<N1; i++)
	    {
	      phi01=phi01+m11array(i)*m11array(i);
	      phi02=phi02+m21array(i)*m21array(i);
	      phi03=phi03+m31array(i)*m31array(i);
	    }
	  //printf("%lf %lf %lf\n",phi01,phi02,phi03);
	  for(int m=1; m<N1; m++)
	    {
	      m11acf(m)=0.0;
	      m21acf(m)=0.0;
	      m31acf(m)=0.0;
	      for(int i=0; i<N1; i++)
		{
		  m11acf(m)=m11acf(m)+m11array(i+m)*m11array(i);
		  m21acf(m)=m21acf(m)+m21array(i+m)*m21array(i);
		  m31acf(m)=m31acf(m)+m31array(i+m)*m31array(i);
		}
	      //printf("%lf %lf %lf %lf\n",(double)m*dt,(m11acf(m-1))/(phi01),(m21acf(m-1))/(phi02),(m31acf(m-1))/(phi03));
	      if(m11shuki(samplecounter)==0.0 && m*dt>=10.0)
		{
		  if((m11acf(m-1))/(phi01)>=0.7)
		    {
		      if((m11acf(m-2))/(phi01)<(m11acf(m-1))/(phi01) && (m11acf(m-1))/(phi01)>(m11acf(m))/(phi01))
			{
			  m11shuki(samplecounter)=(m-1)*dt;
			}
		    }
		}
	      if(m21shuki(samplecounter)==0.0 && m*dt>=10.0)
		{
		  if((m21acf(m-1))/(phi02)>=0.7)
		    {
		      if((m21acf(m-2))/(phi02)<(m21acf(m-1))/(phi02) && (m21acf(m-1))/(phi02)>(m21acf(m))/(phi02))
			{
			  m21shuki(samplecounter)=(m-1)*dt;
			}
		    }
		}
	      if(m31shuki(samplecounter)==0.0 && m*dt>=10.0)
		{
		  if((m31acf(m-1))/(phi03)>=0.7)
		    {
		      if((m31acf(m-2))/(phi03)<(m31acf(m-1))/(phi03) && (m31acf(m-1))/(phi03)>(m31acf(m))/(phi03))
			{
			  m31shuki(samplecounter)=(m-1)*dt;
			}
		    }
		}
	    }
	  //printf("%lf %lf %lf %lf\n",km,m11shuki(0),m21shuki(0),m31shuki(0));
	} //for samplecounter
      
      double average1=0.0;
      double average2=0.0;
      double average3=0.0;
      /*
      double variance1=0.0;
      double variance2=0.0;
      double variance3=0.0;
      */
      for(int samplecounter=0; samplecounter<sample; samplecounter++)
	{
	  average1=average1+m11shuki(samplecounter);
	  average2=average2+m21shuki(samplecounter);
	  average3=average3+m31shuki(samplecounter);
	}
      average1=average1/((double)sample);
      average2=average2/((double)sample);
      average3=average3/((double)sample);
      /*
      for(int samplecounter=0; samplecounter<sample; samplecounter++)
	{
	  variance1=variance1+(m11shuki(samplecounter)-average1)*(m11shuki(samplecounter)-average1);
	  variance2=variance2+(m21shuki(samplecounter)-average2)*(m21shuki(samplecounter)-average2);
	  variance3=variance3+(m31shuki(samplecounter)-average3)*(m31shuki(samplecounter)-average3);
	}
      variance1=variance1/((double)sample);
      variance2=variance2/((double)sample);
      variance3=variance3/((double)sample);
      */
      //printf("%lf %lf %lf %lf %lf %lf %lf\n",km,average1,sqrt(variance1),average2,sqrt(variance2),average3,sqrt(variance3));
      double sa=(average1-average2)*(average1-average2)+(average2-average3)*(average2-average3)+(average3-average1)*(average3-average1);
      printf("%lf %lf %lf %lf %lf %lf\n",km,alpha,sa,average1,average2,average3);

    } //for kmcounter
  return 0;
} //int main()
