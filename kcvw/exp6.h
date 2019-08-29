const char * prefix = "change_in_mean_non_gauss";

const float mu1 = 0;
const float sigma1 = 3; // <--- The std will be 3 and the variance will be 9

const float p_change = 0.02;

const float mu2 = 1;
const float sigma2 = 3; // <-- The std will be 3 and the variance will be 9

const int nref =  10; 

const int famax = 1000;

const float h_delta = 1.;
const int w1_h0 = -10;
const int w5_h0 = -10 ;
const int w10_h0 = -10 ;
const int w15_h0 = -10 ;
const int w20_h0 = -10 ;
const int w40_h0 = -10;

const double cu_delta = 2;

double k(double x, double y, double alpha=2, double sigma2=2){
 return pow( fabs(x) , alpha ) + pow( fabs(y) , alpha) - pow( fabs(x-y) , alpha);
}
