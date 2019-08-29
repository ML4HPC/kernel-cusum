#include <stdio.h>
#include <random>
#include <assert.h>
#include <vector>
#include <list>

#include <string>
#include <ctime>
#include <omp.h>
#include <stdarg.h>

using namespace std;

/* 
  Comparison of ~~kernel~~ CUSUM and window methods.
  Thomas Flynn. 6/01/2018

  Compile with Makefile

  Run with (for example with 2 threads)
  $ export OMP_NUM_THREADS=2
  $ ./kcvw
 */

//Specify which experiment to run by uncommenting the appropriate line:
#ifdef EXP1
#include "exp1.h"
#endif

#ifdef EXP2
#include "exp2.h"
#endif

#ifdef EXP3
#include "exp3.h"
#endif

#ifdef EXP4
#include "exp4.h"
#endif

#ifdef EXP5
#include "exp5.h"
#endif

#ifdef EXP6
#include "exp6.h"
#endif

#ifdef EXP7
#include "exp7.h"
#endif

#ifdef EXP8
#include "exp8.h"
#endif

default_random_engine generator;

const int None = -1;

double h(double x1, double y1, double x2, double y2){
  return k(x1,x2) + k(y1,y2) - k(x1,y2) - k(x2,y1);
}

class normal_variate{
public:
  
  normal_variate():n(0){}
  
  float sample(){
    
    return distribution(generator);
  }
 
  normal_distribution<double> distribution;
  int n;
};

class reference{
public:
  float sample(){
    return sigma1*nv.sample() + mu1;
  }
 
  normal_variate nv;
};

class sequence_generator{
public:
  sequence_generator(int _change_point=None):change_point(_change_point){
    t=0;
  }
  
  float next(){
    float ret;
    
    if ( (t < change_point) or (change_point == None)){
      //printf("Pre-change\n");
      ret = ref.sample();
    }
    else{
      ret = sigma2*nv.sample() + mu2;
    }
    t++;
    return ret;
  }

  int change_point;
  reference ref;
  normal_variate nv;
  int t;
};

template<int width>
class window {
public:
  window():nref(40){
    for(int i=0;i<width;i++) history[i]=0;
    val = 0;
    xprev = ref.sample(); //How to initialize?
  }

  float step(float x){
    
    int i;

    val = 0;
    for(i=0;i<width-1;i++){
      history[i] = history[i+1];
      val += history[i];
    }
    
    history[i] = 0;
    
    for(int n=0;n<nref;n++)
      history[i] += h(xprev, ref.sample(), x, ref.sample());
    
    history[i] /= nref;
    
    xprev = x;
    val += history[i];
    
    if( isnan(val) ){
      printf("wtf x was %f\n",x);
      assert(not isnan(val));
    }
    return val;
  }
  
  float xprev;
  float val;
  float history[width];
  int nref;
  reference ref;
};

class cusum {
public:
  cusum():nref(40){
     val = 0;
     xprev = ref.sample();
  }

  float step(float x){

    float hval = 0;
    
    for(int i=0;i<nref;i++)
      hval += h(xprev, ref.sample(), x, ref.sample());
	
    val = max(0., val + hval/nref - cu_delta );
    xprev = x;
    return val;
  }
  int nref;
  float xprev;
  float val;
  reference ref;
};

template<int nref=1>
class ksprt_lo {
public:
  ksprt_lo():val(0),m(0),kxx(0),kyy(0),kxy(0){}
  
  float step(float x){
    
    for(int i=0;i<nref;i++)
      refprev[m*nref + i] = ref.sample();

    for(int i=0; i<m; i++){
      
      kxx += 2.0*k(x,xprev[i]); //new x with prev x

      for(int j=0;j<nref;j++)
	for(int l=0;l<nref;l++) //new y with prev y
	  kyy += 2.0*k(refprev[m*nref + j] ,
		       refprev[i*nref + l]);
      	

      for(int l=0;l<nref;l++)
	kxy += k(x,refprev[i*nref + l]); //x w prev y
      
      for(int j=0;j<nref;j++)
	kxy += k(refprev[m*nref + j],xprev[i]); 
      
    }

    for(int j=0;j<nref;j++){
      for(int l=0;l<nref;l++){
	if (l == j) continue;
	kyy += k(refprev[m*nref + j] , //new y with previous y
		     refprev[m*nref + l]);
      }
    }

    for(int i=0;i<nref;i++)
      kxy += k(x,refprev[m*nref + i ]); //new x with new y
    
    xprev[m] = x;
  
    m++;
    
    double coeff = pow(m,1.);

    if (nref*m == 1)
      val = 0;
    else
      val = 
	(coeff)*kxx/(m*m) + 
	(coeff)*kyy/(nref*m*(nref*m-1)) - 
	2.0*(coeff)*kxy/(m*m*nref); //we are multiplying by m-1.

    //printf("Val is %f\n",val);
    if(m > 900000){
      printf("OVER 900000!\n");
    }
    
    if (val < 0){
      val = 0 ;
      m = 0;
      kxx = 0;
      kyy = 0;
      kxy = 0;
    }
    return val;
  }
  float refprev[1000000];
  float xprev[1000000];
  
  float val;
  
  float kxx;
  float kyy;
  float kxy;
  
  int m;
  reference ref;
};


template<int nref=1>
class ksprt {  
public:
  ksprt():val(0),m(0){}

  float step(float x){

    for(int i=0;i<nref;i++){
      refprev[m*nref + i] = ref.sample();
    }
    
    for(int i=0; i<m; i++){
      val += k(x,xprev[i]);
      
      for(int j=0;j<nref;j++)
	for(int l=0;l<nref;l++)
	  val += k( refprev[m*nref + j], 
		    refprev[i*nref + l])/(nref*nref);
	

      for(int l=0;l<nref;l++)
	val -= k(x,refprev[i*nref + l])/nref;

      for(int j=0;j<nref;j++)
	val -= k(refprev[m*nref + j], xprev[i])/nref;
    }
    //printf("Val %f\n",val);
    
    xprev[m] = x;
    
    m++;
    if(m > 900000){
      printf("OVER 900000!\n");
    }
    if (val < 0){
      val = 0 ;
      m = 0;
    }

    if (m > 0)
      return val / m;
    //printf("yto\n");
    return 0;
  }
  float refprev[1000000];
  float xprev[1000000];
  
  float val;
  int m;
  reference ref;
};

class trivial {
public:
  trivial(){
     val = 0;
  }

  float step(float){
    
    val += 0.01;
  
    return val;
  }


  float val;

};


class average {
public:
  average(){
    cumulative = 0;
    n = 0;
  }
  
  void tally(float val){
    cumulative += val;
    n++;
  }

  float avg(){
    if (n == 0)
      return 0;
    else
      return cumulative/n;
  }
  
  float cumulative;
  int n;
  
};


class profile {
public:
  profile(){
    prev = time(0);
  }

  void step(){
    time_t cur = time(0);
    //difftime returns seconds.
    //multiplying by 1k yields milliseconds.
    times.tally( difftime( cur, prev) * 1000.0 );
    prev = cur;
  }

  float avg(){ 
    return times.avg();
  }
  time_t prev;
  average times;
};

class csv {
public:
  csv(const char * fname = NULL, const char * prefix = "."):f(NULL){
    
    if(fname){
      char totalname[255];
      sprintf(totalname,"%s/%s",prefix,fname);
      f = fopen(totalname,"w");
      printf("Opened the file %s\n",totalname);
      assert(f);
    }
    else{
      printf("not opening the file...\n");
    }
  }

  void row(const vector<float> & data){
    if (not f) return;
    if (data.size() == 0) return;

    int i;
    for(i = 0; i < data.size()-1;i++)
      fprintf(f, "%f,",data[i]);
    fprintf(f, "%f\n",data[i]);

    fflush(f);
  }

  ~csv(){
    if (f)
      fclose(f);
  }
  FILE * f;
};

struct cd_props {
  float ttfa ;
  float delays;
};

template<class Statistic>
int detect(float h, int cpt){
  sequence_generator seq(cpt);
  Statistic stat = Statistic();
  int t = 0;
  while (stat.step( seq.next() ) < h) t++;
  
  //printf("Final stat was %f\n",stat.val);
  return t; 
}

template<class Statistic>
float false_alarms(float h, int n ){
  
  //false alarm
  average false_alarms;
  profile prof;

  for(int i=0; i<n ;i++){
    int t = detect<Statistic>  (h, None);
    false_alarms.tally(t);
    prof.step();
    //printf("Got FA=%d at i=%d",t,i);
  }
  printf("false alarms took %f ms on average \n",prof.avg());
  return false_alarms.avg();
}

template<class Statistic>
float delays(float h, int n){
 
  
  average delays;
  geometric_distribution<int> geom(p_change);
  for(int i=0; i < n; i++){
    int cpt =  geom(generator);
    //cpt = 50;
    int t = detect< Statistic > (h, cpt);
    //printf("The cpt was %d and the t was %d and the h was %f\n", cpt,t,h);
    if (t > cpt) delays.tally(t-cpt);
  }

  return delays.avg();
}

template<class Statistic>
cd_props false_alarms_and_delays(float h, int n){

  
  float fas = false_alarms<Statistic>(h, n);
 
  float ds = delays<Statistic>(h, n);
  
  return { fas, ds };
}

char *  str(const char* format, ...)
{
  char * chr_array = new char[255];
  va_list argptr;
  va_start(argptr, format);
  vsprintf(chr_array, format, argptr);
  va_end(argptr);
  return chr_array;
}

// Find the relationship between false alarms and delays by running through many threshold values
template<class Statistic>
void calculate_response(const char * fname = NULL, 
			float h0=0, float h_delta=1, int famax=200, int seed=None, int n = 3000){
  if (seed != None)
    generator.seed(seed);

  printf("Doing response %s\n",(fname)?fname:"(no file)" );
  //delay
  cd_props res; res.ttfa = res.delays = 0;
  float h = h0;
  
  csv outfile(fname,prefix);
  
  while (res.ttfa < famax){
    
    profile prof;
    
    res = false_alarms_and_delays< Statistic >(h,n);
    prof.step();

    printf("at h = %f False alarm average: %f delay average: %f time was %f \n", 
	   h, res.ttfa, res.delays,prof.avg());

    outfile.row({h,res.ttfa,res.delays});
    h += h_delta;
  }
  printf("*************** Done with %s\n",fname);
}

// Find the threshold h resulting in a desired false-alarm rate.
template<class Statistic>
void solve_for_ttfa(float h0=0, float hmax = 100, int fatgt=500, int seed=None){
  
  if (seed != None) generator.seed(seed);
  
  float fa0 = false_alarms<Statistic>(h0);
  float famax = false_alarms<Statistic>(hmax);
 
  if ((fa0 > fatgt) or (famax < fatgt)){
    printf("Invalid initial conditions: h0 %f yields fa %f and hmax %f yields fa %f",
	   h0,fa0,hmax,famax);
    return;
  }

  printf("Interval is [%f, %f]. ", h0,hmax);
  float hmid = (h0 + hmax)/2;  
  float famid = false_alarms<Statistic>(hmid);
  printf("midpt h = %f yields false alarm average %f.\n", hmid, famid);

  if (famid < fatgt)
    h0 = hmid;
  else
    hmax = hmid;

 
  while (true){
   
    printf("Interval is [%f, %f]\t", h0,hmax);
    hmid = (h0 + hmax)/2;

    famid = false_alarms<Statistic>(hmid);

    if (famid < fatgt)
      h0 = hmid;
    else
      hmax = hmid;
  
    printf("at h = %f False alarm average: %f\n", hmid, famid);

  }

}

template<class Statistic>
void kex(int h){
  geometric_distribution<int> geom(0.02);
  int cpt =  500; //geom(generator);
  printf("Change point is %d\n",cpt);
  sequence_generator seq(cpt);
  Statistic stat = Statistic();
  csv seqcsv("k_seq.csv");
  csv statcsv("k_stat.csv");
  
  float seqval =  seq.next();
  float statval = stat.step(seqval);
  int t=0;
  while (t < 1000 ){
    seqcsv.row({seqval});
    statcsv.row({statval});

    seqval =  seq.next();
    statval = stat.step(seqval);
  
    t++;
  }
}

void error(const char * msg){
  printf("%s\n",msg);
  exit(1);
}

int main(int argc, char ** argv){
  int seed = time(NULL);
  srand48(seed);
  srand(seed);
  
  //solve_for_ttfa<window<20> > (-4, 5, 500); //Result is around 7.278
  //solve_for_ttfa<cusum > (-4, 10, 500); //Result is around 9.254

  //Calculate delays at TTFA = 500.
  //float dsw = delays<window<20> >(7.278);
  //float dsc = delays<cusum>(9.254);
  //printf("Window delays: %f CUSUM delays %f\n", dsw,dsc);
  
  //calculate_response< trivial >("t.csv", 0, 1000,seed);

  /*  cd_props res;
  float h = 2;
  res = false_alarms_and_delays< cusum >(h);
  printf("at h = %f False alarm average: %f delay average: %f\n",h , res.ttfa, res.delays);
  return 0;
  */

  //  kex<ksprt>(3);
  // return 0;
  if(argc < 2){
    error("no thread argument specified.");
  }
  
  int tid = atoi(argv[1]);

  //  calculate_response< ksprt<10> >("test.csv",0.8, 0.2 , 200);
  //return 0;

  printf("Tid is %d\n",tid);
  generator.seed(time(NULL) + tid);


  //calculate_response< ksprt_lo<10> >(str("k10_lo_new_%d.csv", tid),0, 0.025 , famax);
  //return 0;

  //  calculate_response< ksprt_lo<1> >(str("k5_lo_%d.csv", tid),0, 0.1 , famax);
  //return 0;

  //calculate_response< ksprt_lo<5> >(str("k5_lo_%d.csv", tid),0, 0.05 , famax);
  //return 0;

  if(tid % 3 == 0)//default h_delta is 0.4
     calculate_response< ksprt<10>  >(str("k10s_%d.csv", tid),0, 0.025, famax);
  else if (tid % 3 == 1)
    calculate_response< ksprt_lo<10> >(str("klo10s_%d.csv", tid),0, 0.025 , famax);
  else if(tid % 3 == 2)
    calculate_response< window<10>  >(str("w10s_%d.csv", tid) ,0, 0.025, famax);
  /*else if(tid % 6 == 3)
    ;    calculate_response< ksprt<1>  >(str("k1hi_%d.csv", tid),0, 0.1, famax);
  //calculate_response< window<40> >("w40.csv", w10_h0 , famax);
  else if(tid % 6 == 4)
    calculate_response< ksprt<5>  >(str("k5hi_%d.csv", tid),0, 0.1, famax);
  //calculate_response< window<15> >("w15.csv", w15_h0, famax);
  else if(tid % 6 == 5)
    calculate_response< ksprt<10>  >(str("k10hi_%d.csv", tid),0, 0.1, famax);
     */

  //calculate_response< window<20> >("w20.csv",w20_h0, famax);
  else if(tid == 6)
    ;//calculate_response< window<40> >("w40.csv",w40_h0, famax );
    

  return 0;
}
