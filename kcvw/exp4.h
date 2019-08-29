//Differs from exp3 by more reference data
const char * prefix = "small_change_in_std_more_ref";
const float mu1 = 0;
const float sigma1 = 1;  // <--- The std will be 1 and the variance will be 1.
const float p_change = 0.02;
const float mu2 = 0 ;
const float sigma2 = 1.1; // <-- The std will be 1.5 and the variance will be 1.5^2
const int nref = 20; 

const int famax = 1000;

const float h_delta = 1;
const int w1_h0 = -2;
const int w5_h0 = -2;
const int w10_h0 = -2;
const int w15_h0 = -2 ;
const int w20_h0 = -2;
const int w40_h0 = -2;

const double cu_delta = 2;

double k(double x, double y, double alpha=2, double sigma2=2){
  return pow( fabs(x) , alpha ) + pow( fabs(y) , alpha) - pow( fabs(x-y) , alpha);
}
