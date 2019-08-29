#include <stdio.h>
#include <random>
#include <assert.h>
#include <vector>
#include <string>
#include <time.h>

/* 
  Comparison of  CUSUM and window methods.
  Thomas Flynn. 5/31/2018


  Compile with g++ -std=c++11 cusum-vs-window.cpp -o cvw
  Run with ./cvw
 */

using namespace std;

default_random_engine generator;

double ratio(float x,float mu1 = 0, float mu2=1, float sigma=3){
  return ((mu2 - mu1)/(2*sigma))*(2*x - (mu1 + mu2));
}
  
template<int width>
class window {
public:
  window(){
    for(int i=0;i<width;i++) history[i]=0;
    val = 0;
  }

  float step(float x){
    int i;
    val = 0;
    for(i=0;i<width-1;i++){
      history[i] = history[i+1];
      val += history[i];
    }
    history[i] = ratio(x);
    val += history[i];
    return val;
  }
  
  float val;
  float history[width];
};


class cusum {
public:
  cusum(){
     val = 0;
  }

  float step(float x){
    val = max(0., val + ratio(x));
    
    return val;
  }
  
  float val;
};

const int None = -1;

class normal_variate{
public:
  float sample(){
    return distribution(generator);
  }
 
  normal_distribution<double> distribution;
};

template<int change_point>
class sequence_generator{
public:
  sequence_generator(){
    
    t=0;
  }
  
  float next(){
    float ret;
    
    if ( (t < change_point) or (change_point == None)){
      //printf("Pre-change\n");
      ret = nv.sample()*3;
    }
    else{
      ret = nv.sample()*3 + 1;
    }
    t++;
    return ret;
  }
  
  normal_variate nv;
  int t;
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

class csv {
public:
  csv(const char * fname){
    f = fopen(fname,"w");
    assert(f);
  }

  void row(const vector<float> & data){
    int i;
    for(i = 0; i<data.size()-1;i++)
      fprintf(f, "%f,",data[i]);
    fprintf(f, "%f\n",data[i]);
  }

  ~csv(){
    fclose(f);
  }
  FILE * f;
};

template<class Sequence, class Statistic>
int detect(float h){
  Sequence seq = Sequence();
  Statistic stat = Statistic();
  int t = 0;
  while (stat.step( seq.next() ) < h) t++;
  return t; 
}

struct cd_props {
  float ttfa ;
  float delays;
};


template<class Statistic>
float false_alarms(float h, int n = 30000){
  
  //false alarm
  average false_alarms;
  for(int i=0; i<n ;i++){
    //printf("Another one...\n");
    int t = detect< sequence_generator<None>, Statistic>  (h);
    false_alarms.tally(t);
  }
 
  return false_alarms.avg();
}

template<class Statistic>
float delays(float h, int n = 30000){
  const int cpt=50;
  
  average delays;
  
  for(int i=0; i < n; i++){

    int t = detect<sequence_generator<cpt>, Statistic > (h);
    
    if (t > cpt) delays.tally(t-cpt);
  }

  return delays.avg();
}

template<class Statistic>
cd_props false_alarms_and_delays(float h, int n = 30000){
  
  float fas = false_alarms<Statistic>(h, 30000);
  float ds = delays<Statistic>(h,30000);
  
  return { fas, ds };
}

string str(const char * s){
  return string(s);
}

string str(int w){
  char buff[255];
  sprintf(buff,"%d", w);
  return string(buff);
}

// Find the relationship between false alarms and delays
template<class Statistic>
void calculate_response(const char * fname,float h0=0,int famax=200, int seed=None){
  if (seed != None)
    generator.seed(seed);

  normal_variate nv;
  printf("Some normal samples: %f, %f, %f", nv.sample(),nv.sample(),nv.sample());
  
  printf("Doing response %s\n",fname);
  //delay
  cd_props res; res.ttfa = res.delays = 0;
  float h = h0;
  
  csv outfile(fname);
  while (res.ttfa < famax){
    res = false_alarms_and_delays< Statistic >(h);
    printf("at h = %f False alarm average: %f delay average: %f\n", h, res.ttfa, res.delays);
    h += 0.1;
    outfile.row({h,res.ttfa,res.delays});
  }

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


int main(){
  int seed = time(NULL);

  //solve_for_ttfa<window<20> > (-4, 10, 500); //Result is around 7.278
  //solve_for_ttfa<cusum > (-4, 10, 500); //Result is around 9.254
  //float dsw = delays<window<20> >(7.278);
  //float dsc = delays<cusum>(9.254);
  //printf("Window delays: %f CUSUM delays %f\n", dsw,dsc);
  calculate_response< window<1> >("w1.csv", -1, 1000,seed);
  calculate_response< window<5> >("w5.csv", -1, 1000);
  calculate_response< window<10> >("w10.csv", -3, 1000);
  calculate_response< window<15> >("w15.csv", -3, 1000);
  calculate_response< window<20> >("w20.csv",-4, 1000);
  calculate_response< window<40> >("w40.csv",-4, 1000, seed );
  calculate_response< cusum >("c.csv",0, 1000, seed);

  return 0;
}
