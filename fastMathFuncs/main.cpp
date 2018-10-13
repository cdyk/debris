#include <cmath>
#include <iostream>

#include <emmintrin.h>  // SSE2

float asin_ss(float _x)
{
  static const float signBit = -0.f;
  __m128 sign = _mm_load_ss(&signBit);

  __m128 x = _mm_set_ss(_x);
  __m128 abs_x = _mm_andnot_ps(sign, x);
 
#if 0
  // Max error < 2e-5
  static const float C3 = -0.0187293f;
  static const float C2 =  0.0742610f;
  static const float C1 = -0.2121144f;
  static const float C0 =  1.5707288f;
  __m128 r = _mm_mul_ss(_mm_load_ss(&C3), abs_x);
  r = _mm_add_ss(r, _mm_mul_ss(_mm_load_ss(&C2), abs_x));
  r = _mm_add_ss(r, _mm_mul_ss(_mm_load_ss(&C1), abs_x));
  r = _mm_add_ss(r, _mm_load_ss(&C0));
#elif 1
  // Max error < 2e-8
  static const float C7 = -0.0012624911f;
  static const float C6 =  0.0066700901f;
  static const float C5 = -0.0170881256f;
  static const float C4 =  0.0308918810f;
  static const float C3 = -0.0501743046f;
  static const float C2 =  0.0889789874f;
  static const float C1 = -0.2145988016f;
  static const float C0 =  1.5707963050f;
  __m128 r = _mm_mul_ss(_mm_load_ss(&C7), abs_x);
  r = _mm_add_ss(r, _mm_mul_ss(_mm_load_ss(&C6), abs_x));
  r = _mm_add_ss(r, _mm_mul_ss(_mm_load_ss(&C5), abs_x));
  r = _mm_add_ss(r, _mm_mul_ss(_mm_load_ss(&C4), abs_x));
  r = _mm_add_ss(r, _mm_mul_ss(_mm_load_ss(&C3), abs_x));
  r = _mm_add_ss(r, _mm_mul_ss(_mm_load_ss(&C2), abs_x));
  r = _mm_add_ss(r, _mm_mul_ss(_mm_load_ss(&C1), abs_x));
  r = _mm_add_ss(r, _mm_load_ss(&C0));
#endif

  static const float one = 1.f;
  __m128 q = _mm_sub_ss(_mm_load_ss(&one), abs_x);
#if 0
  // lat/thru 25/16
  q = _mm_sqrt_ss(q);
#else
  // Replace high-precision sqrt with 1/(1/sqr) combo
  // lat/thru -> 3/3 + 3/3
  q = _mm_rcp_ss(_mm_rsqrt_ss(q));
#endif

  const float piTwo = 1.5707963267948966f;
  r = _mm_sub_ss(_mm_load_ss(&piTwo), _mm_mul_ss(q, r));

  // copy sign from x
  r = _mm_or_ps(r, _mm_and_ps(x, sign));

  float rv;
  _mm_store_ss(&rv, r);
  return rv;
}



#include <cmath>
#include <iostream>

#include <emmintrin.h>  // SSE2

float atanApprox(float x)
{
  float x_x = x*x;
  float r = 0.0208351f*x_x;
  r = (r - 0.0851330f)*x_x;
  r = (r + 0.1801410f)*x_x;
  r = (r - 0.3302995f)*x_x;
  r = (r + 0.9998660f)*x;
  return r;
}


float atan2Approx(float y, float x)
{
  const float piTwo = 1.5707963267948966f;
  const float pi = 3.1415926535897931f;
  
  bool pq = std::abs(y) < std::abs(x);
  
  float num = pq ? y : x;
  float den = pq ? x : y;
  float t0 = atanApprox(num/den);
  float atan = pq ? t0 : -t0;

  bool px = 0 <= x;
  bool py = 0 <= y;
  float t2 = px ? 0.f : pi;
  float t3 = pq ? t2 : piTwo;
  float shift = py ? t3 : -t3;

  return atan + shift;
}

float atan2_ss(float _y, float _x)
{
  static const float piTwo = 1.5707963267948966f;
  static const float pi = 3.1415926535897931f;
  static const float C9 =  0.0208351f;
  static const float C7 = -0.0851330f;
  static const float C5 =  0.1801410f;
  static const float C3 = -0.3302995f;
  static const float C1 =  0.9998660f;
  static const float signBit = -0.f;

  __m128 sign = _mm_load_ss(&signBit);
  __m128 x = _mm_set_ss(_x);
  __m128 y = _mm_set_ss(_y);
  __m128 abs_x = _mm_andnot_ps(sign, x);
  __m128 abs_y = _mm_andnot_ps(sign, y);
  __m128 pq = _mm_cmplt_ss(abs_y, abs_x);

  __m128 num = _mm_or_ps( _mm_and_ps(pq, y), _mm_andnot_ps(pq, x));
  __m128 den = _mm_or_ps( _mm_and_ps(pq, x), _mm_andnot_ps(pq, y));
  __m128 t = _mm_div_ss(num, den);

  __m128 t_t = _mm_mul_ss(t, t);

  __m128 r = _mm_mul_ss(_mm_set_ss(C9), t_t);
  r = _mm_mul_ss(_mm_add_ss(r, _mm_load_ss(&C7)), t_t);
  r = _mm_mul_ss(_mm_add_ss(r, _mm_load_ss(&C5)), t_t);
  r = _mm_mul_ss(_mm_add_ss(r, _mm_load_ss(&C3)), t_t);
  r = _mm_mul_ss(_mm_add_ss(r, _mm_load_ss(&C1)), t);
  
  r = _mm_xor_ps(r, _mm_andnot_ps(pq, sign));
  
  __m128 t2 = _mm_andnot_ps(_mm_cmple_ss(_mm_setzero_ps(), x), _mm_set_ss(pi));     // t2    = 0 <= x ? 0  : pi
  __m128 t3 = _mm_or_ps(_mm_and_ps(pq, t2), _mm_andnot_ps(pq, _mm_set_ss(piTwo)));  // t3    =   pq   ? t2 : pi/2
  __m128 shift = _mm_xor_ps(t3, _mm_and_ps(y, sign));                               // shift = 0 <= y ? t3 : - t3
 

  float rv;
  _mm_store_ss(&rv, _mm_add_ss(r, shift));
  
  return rv;
}



int main()
{

  for (int i = 0; i < 100; i++) {
    float x = (2.f*i) / 100 - 1.f;
    std::cout << x << '\t' << asin_ss(x) << "\t" << asin(x) << "\n";
  }

  for(int j=0; j<100; j++) {
    float y = (2.f*j)/100 - 1.f;

      for(int i=0; i<100; i++) {
  float x = (2.f*i)/100 - 1.f;

  std::cout << x << '\t' << y << '\t' << atan2_ss(y,x) << "\t" << atan2(y,x) << "\n";
      }
      std::cout << "\n";
  }
  return 0;
}
