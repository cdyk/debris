#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>

#include "xmmintrin.h"
#include "immintrin.h"

typedef float Vec3[3];

void vecSub(Vec3 & D, const Vec3& A, const Vec3& B)
{
  D[0] = A[0] - B[0];
  D[1] = A[1] - B[1];
  D[2] = A[2] - B[2];
}

void vecCross(Vec3& D, const Vec3& A, const Vec3& B)
{
  D[0] = A[1] * B[2] - A[2] * B[1];
  D[1] = A[2] * B[0] - A[0] * B[2];
  D[2] = A[0] * B[1] - A[1] * B[0];
}

float vecDot(const Vec3& A, const Vec3& B)
{
  return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}

float vecDist(const Vec3& A, const Vec3& B)
{
  float t[3];
  vecSub(t, A, B);
  return vecDot(t, t);
}

template<typename T>
T mask(bool cond, const T V)
{
  return cond ? V : T(0);
}

//N[0] = u[1]v[2] - u[2]v[1];
//N[1] = u[2]v[0] - u[0]v[2];
//N[2] = u[0]v[1] - u[1]v[0];

unsigned nearestPointOnTriangle(float(&Q)[3], const float(&A)[3], const float(&B)[3], const float(&C)[3], const float(&P)[3])
{
  // Algorithm:
  // - For each edge, form an outwards plane NxAB, NxBC, and NxCA.
  // - Distance from these planes, e.g. <AP,NxAB> forms barycentric coords when they are normalized.
  // - These half-planes also forms the boundary of the Voronoi-cells of the vertices and edges of the triangle.
  // - By inspecting the sign of the plane distances, we can determine which region P is in:
  // 
  //        \110=6 |
  //         \CA<0 |
  //          \BC<0|
  //           \   |
  //            \  |
  //             \ |
  //              \|
  //               C
  //          CA<0 |\   BC<0
  //         100=4 | \  010=2
  //               |  \
  //        -------A---B------------
  //          AB<0 |AB<0\ AB<0
  //          CA<0 |     \ BC<0
  //         101=5 | 001=1\ 011=3
  //
  // - If all distances are positive, the closest point is inside the triangle,
  //   and we use the barycentric coords to find the position.
  // - If one distance is positive and two negative, the closest point is a
  //   corner, and we use that position.
  // - If two distances are positive and one negative, the closest point is on
  //   an edge. If on edge AB, the ratio of <AP,AB> and <PB,AB> (P orthognoally
  //   projected onto AB) gives the position.

  float AP[3];
  float BP[3];
  float CP[3];
  float N[3];
  float AB[3];
  float BC[3];
  float CA[3];
  for (unsigned i = 0; i < 3; i++) {
    AP[i] = P[i] - A[i];
    BP[i] = P[i] - B[i];
    CP[i] = P[i] - C[i];

    AB[i] = B[i] - A[i];
    BC[i] = C[i] - B[i];
    CA[i] = A[i] - C[i];
  }

  // Calc barycentric coords.
  N[0] = AB[1] * BC[2] - AB[2] * BC[1];
  N[1] = AB[2] * BC[0] - AB[0] * BC[2];
  N[2] = AB[0] * BC[1] - AB[1] * BC[0];

  float NxAB[3];
  NxAB[0] = N[1] * AB[2] - N[2] * AB[1];
  NxAB[1] = N[2] * AB[0] - N[0] * AB[2];
  NxAB[2] = N[0] * AB[1] - N[1] * AB[0];
  float AP_dot_NxAB =
    AP[0] * NxAB[0] +
    AP[1] * NxAB[1] +
    AP[2] * NxAB[2];

  float NxBC[3];
  NxBC[0] = N[1] * BC[2] - N[2] * BC[1];
  NxBC[1] = N[2] * BC[0] - N[0] * BC[2];
  NxBC[2] = N[0] * BC[1] - N[1] * BC[0];
  float BP_dot_NxBC =
    BP[0] * NxBC[0] +
    BP[1] * NxBC[1] +
    BP[2] * NxBC[2];

  float NxCA[3];
  NxCA[0] = N[1] * CA[2] - N[2] * CA[1];
  NxCA[1] = N[2] * CA[0] - N[0] * CA[2];
  NxCA[2] = N[0] * CA[1] - N[1] * CA[0];
  float CP_dot_NxCA =
    CP[0] * NxCA[0] +
    CP[1] * NxCA[1] +
    CP[2] * NxCA[2];


  uint32_t region = (CP_dot_NxCA < 0.f ? 4 : 0) | (BP_dot_NxBC < 0.f ? 2 : 0) | (AP_dot_NxAB < 0.f ? 1 : 0);
  float w_a = 0.f;
  float w_b = 0.f;
  float w_c = 0.f;

  float* e0 = nullptr;
  float* e1 = nullptr;
  float* e2 = nullptr;

  switch (region)
  {
  case 0:   // 000 - inside ab, bc, ca
    w_a = BP_dot_NxBC;
    w_b = CP_dot_NxCA;
    w_c = AP_dot_NxAB;
    break;
   
  case 3:   // 011 - outside ab, bc; inside ca;
    w_b = 1.f;
    break;
  case 5:   // 101 - outside ca, ab; inside bc;
    w_a = 1.f;
    break;
  case 6:   // 110 - outside ca, bc; inside ab
    w_c = 1.f;
    break;

  case 1:   // 001 - outside ab; inside bc, ca;
    w_a = -(BP[0] * AB[0] + BP[1] * AB[1] + BP[2] * AB[2]);
    w_b = AP[0] * AB[0] + AP[1] * AB[1] + AP[2] * AB[2];
    break;
  case 2:   // 010 - outside bc; inside ab, ca;
    w_b = -(CP[0] * BC[0] + CP[1] * BC[1] + CP[2] * BC[2]);
    w_c = BP[0] * BC[0] + BP[1] * BC[1] + BP[2] * BC[2];
    break;
  case 4:   // 100 - outside ca; inside bc, ab;
    e0 = AP;
    w_c = -(AP[0] * CA[0] + AP[1] * CA[1] + AP[2] * CA[2]);
    w_a = CP[0] * CA[0] + CP[1] * CA[1] + CP[2] * CA[2];
    break;

  case 7:   // 111 - outside ca, bc, ab
    w_a = w_b = w_c = 1.f;
    break;
  default:
    break;
  }

  if (e0) {


  }

  w_a = w_a < 0.f ? 0.f : w_a;
  w_b = w_b < 0.f ? 0.f : w_b;
  w_c = w_c < 0.f ? 0.f : w_c;

  float s = 1.f/(w_a + w_b + w_c);
  for (unsigned i = 0; i < 3; i++) {
    Q[i] = s * (w_a * A[i] + w_b * B[i] + w_c * C[i]);
  }

  return region;
}


unsigned nearestPointOnTriangle2(float(&Q)[3], const float(&A)[3], const float(&B)[3], const float(&C)[3], const float(&P)[3])
{
  float AB[3], BC[3], CA[3];
  vecSub(AB, B, A);
  vecSub(BC, C, B);
  vecSub(CA, A, C);

  float N[3];
  vecCross(N, AB, BC);

  float NxAB[3], NxBC[3], NxCA[3];
  vecCross(NxAB, N, AB);
  vecCross(NxBC, N, BC);
  vecCross(NxCA, N, CA);

  float AP[3], BP[3], CP[3];
  vecSub(AP, P, A);
  vecSub(BP, P, B);
  vecSub(CP, P, C);
  float AP_dot_NxAB = vecDot(AP, NxAB);
  float BP_dot_NxBC = vecDot(BP, NxBC);
  float CP_dot_NxCA = vecDot(CP, NxCA);

  float AP_dot_AB = vecDot(AP, AB);
  float BP_dot_BC = vecDot(BP, BC);
  float CP_dot_CA = vecDot(CP, CA);

  float neg_BP_dot_AB = -vecDot(BP, AB);
  float neg_CP_dot_BC = -vecDot(CP, BC);
  float neg_AP_dot_CA = -vecDot(AP, CA);

  bool AP_dot_AB_p = 0.f <= AP_dot_AB;
  bool BP_dot_BC_p = 0.f <= BP_dot_BC;
  bool CP_dot_CA_p = 0.f <= CP_dot_CA;

  bool AP_dot_NxAB_n = AP_dot_NxAB <= 0.f;
  bool BP_dot_NxBC_n = BP_dot_NxBC <= 0.f;
  bool CP_dot_NxCA_n = CP_dot_NxCA <= 0.f;

  bool neg_AP_dot_CA_n = neg_AP_dot_CA <= 0.f;
  bool neg_BP_dot_AB_n = neg_BP_dot_AB <= 0.f;
  bool neg_CP_dot_BC_n = neg_CP_dot_BC <= 0.f;

  bool edge_ab = AP_dot_AB_p ? AP_dot_NxAB_n && !neg_BP_dot_AB_n : neg_AP_dot_CA_n;
  bool edge_bc = BP_dot_BC_p ? BP_dot_NxBC_n && !neg_CP_dot_BC_n : neg_BP_dot_AB_n; 
  bool edge_ca = CP_dot_CA_p ? CP_dot_NxCA_n && !neg_AP_dot_CA_n : neg_CP_dot_BC_n;

  bool out_ab = edge_ab;
  bool out_bc = !out_ab && edge_bc;
  bool out_ca = !(out_ab || out_bc) && edge_ca;

  bool out_none = !(out_ab || out_bc || out_ca);

  uint32_t region = mask(out_ab, 1) + mask(out_bc, 2) + mask(out_ca, 3);
  float w_a_ = mask(out_ab, neg_BP_dot_AB) + mask(out_ca, CP_dot_CA) + mask(out_none, BP_dot_NxBC);
  float w_b_ = mask(out_ab, AP_dot_AB) + mask(out_bc, neg_CP_dot_BC) + mask(out_none, CP_dot_NxCA);
  float w_c_ = mask(out_bc, BP_dot_BC) + mask(out_ca, neg_AP_dot_CA) + mask(out_none, AP_dot_NxAB);

  float w_a = w_a_ < 0.f ? 0.f : w_a_;
  float w_b = w_b_ < 0.f ? 0.f : w_b_;
  float w_c = w_c_ < 0.f ? 0.f : w_c_;

  float s = 1.f / (w_a + w_b + w_c);
  if (std::isinf(s)) {
    s = w_a = 1.f;
    w_b = w_c = 0.f;
  }

  for (unsigned i = 0; i < 3; i++) {
    Q[i] = s * (w_a * A[i] + w_b * B[i] + w_c * C[i]);
  }

  float pq = vecDist(P, Q);
  float pa = vecDist(P, A);
  float pb = vecDist(P, B);
  float pc = vecDist(P, C);

  assert(pq <= pa * (1 + 1e-5f));
  assert(pq <= pb * (1 + 1e-5f));
  assert(pq <= pc * (1 + 1e-5f));

  return region;
}

unsigned nearestPointOnTriangleSSE(float(&Q_)[3], const float(&A_)[3], const float(&B_)[3], const float(&C_)[3], const float(&P_)[3])
{
  __m128 Ax = _mm_load_ps(A_);  // A.x A.y A.z of triangle 0
  __m128 Ay = _mm_load_ps(A_);  // A.x A.y A.z of triangle 1
  __m128 Az = _mm_load_ps(A_);  // A.x A.y A.z of triangle 2
  __m128 At = _mm_load_ps(A_);  // A.x A.y A.z of triangle 3
  _MM_TRANSPOSE4_PS(Ax, Ay, Az, At);

  __m128 Bx = _mm_load_ps(B_);
  __m128 By = _mm_load_ps(B_);
  __m128 Bz = _mm_load_ps(B_);
  __m128 Bt = _mm_load_ps(B_);
  _MM_TRANSPOSE4_PS(Bx, By, Bz, Bt);

  __m128 Cx = _mm_load_ps(C_);
  __m128 Cy = _mm_load_ps(C_);
  __m128 Cz = _mm_load_ps(C_);
  __m128 Ct = _mm_load_ps(C_);
  _MM_TRANSPOSE4_PS(Cx, Cy, Cz, Ct);

  //vecSub(AB, B, A);
  __m128 ABx = _mm_sub_ps(Bx, Ax);
  __m128 ABy = _mm_sub_ps(By, Ay);
  __m128 ABz = _mm_sub_ps(Bz, Az);

  //vecSub(BC, C, B);
  __m128 BCx = _mm_sub_ps(Cx, Bx);
  __m128 BCy = _mm_sub_ps(Cy, By);
  __m128 BCz = _mm_sub_ps(Cz, Bz);

  //vecSub(CA, A, C);
  __m128 CAx = _mm_sub_ps(Ax, Cx);
  __m128 CAy = _mm_sub_ps(Ay, Cy);
  __m128 CAz = _mm_sub_ps(Az, Cz);

  //vecCross(N, AB, BC);
  __m128 Nx = _mm_fmsub_ps(ABy, BCz, _mm_mul_ps(ABz, BCy));
  __m128 Ny = _mm_fmsub_ps(ABz, BCx, _mm_mul_ps(ABx, BCz));
  __m128 Nz = _mm_fmsub_ps(ABx, BCy, _mm_mul_ps(ABy, BCx));

  // vecCross(NxAB, N, AB);
  __m128 NxABx = _mm_fmsub_ps(Ny, ABz, _mm_mul_ps(Nz, ABy));
  __m128 NxABy = _mm_fmsub_ps(Nz, ABx, _mm_mul_ps(Nx, ABz));
  __m128 NxABz = _mm_fmsub_ps(Nx, ABy, _mm_mul_ps(Ny, ABx));

  //vecCross(NxBC, N, BC);
  __m128 NxBCx = _mm_fmsub_ps(Ny, BCz, _mm_mul_ps(Nz, BCy));
  __m128 NxBCy = _mm_fmsub_ps(Nz, BCx, _mm_mul_ps(Nx, BCz));
  __m128 NxBCz = _mm_fmsub_ps(Nx, BCy, _mm_mul_ps(Ny, BCx));

  //vecCross(NxCA, N, CA);
  __m128 NxCAx = _mm_fmsub_ps(Ny, CAz, _mm_mul_ps(Nz, CAy));
  __m128 NxCAy = _mm_fmsub_ps(Nz, CAx, _mm_mul_ps(Nx, CAz));
  __m128 NxCAz = _mm_fmsub_ps(Nx, CAy, _mm_mul_ps(Ny, CAx));

  // -- per-point loop 

  __m128 signbit = _mm_set_ps1(-0.f);

  __m128 Pt = _mm_loadu_ps(P_);
  __m128 Px = _mm_shuffle_ps(Pt, Pt, _MM_SHUFFLE(0, 0, 0, 0));
  __m128 Py = _mm_shuffle_ps(Pt, Pt, _MM_SHUFFLE(1, 1, 1, 1));
  __m128 Pz = _mm_shuffle_ps(Pt, Pt, _MM_SHUFFLE(2, 2, 2, 2));

  //vecSub(AP, P, A);
  __m128 APx = _mm_sub_ps(Px, Ax);
  __m128 APy = _mm_sub_ps(Py, Ay);
  __m128 APz = _mm_sub_ps(Pz, Az);
  __m128 AP_dot_NxAB = _mm_fmadd_ps(APx, NxABx, _mm_fmadd_ps(APy, NxABy, _mm_mul_ps(APz, NxABz)));
  __m128 AP_dot_AB = _mm_fmadd_ps(APx, ABx, _mm_fmadd_ps(APy, ABy, _mm_mul_ps(APz, ABz)));
  __m128 neg_AP_dot_CA = _mm_xor_ps(signbit, _mm_fmadd_ps(APx, CAx, _mm_fmadd_ps(APy, CAy, _mm_mul_ps(APz, CAz))));

  //vecSub(BP, P, B);
  __m128 BPx = _mm_sub_ps(Px, Bx);
  __m128 BPy = _mm_sub_ps(Py, By);
  __m128 BPz = _mm_sub_ps(Pz, Bz);
  __m128 BP_dot_NxBC = _mm_fmadd_ps(BPx, NxBCx, _mm_fmadd_ps(BPy, NxBCy, _mm_mul_ps(BPz, NxBCz)));
  __m128 BP_dot_BC = _mm_fmadd_ps(BPx, BCx, _mm_fmadd_ps(BPy, BCy, _mm_mul_ps(BPz, BCz)));
  __m128 neg_BP_dot_AB = _mm_xor_ps(signbit, _mm_fmadd_ps(BPx, ABx, _mm_fmadd_ps(BPy, ABy, _mm_mul_ps(BPz, ABz))));

  //vecSub(CP, P, C);
  __m128 CPx = _mm_sub_ps(Px, Cx);
  __m128 CPy = _mm_sub_ps(Py, Cy);
  __m128 CPz = _mm_sub_ps(Pz, Cz);
  __m128 CP_dot_NxCA = _mm_fmadd_ps(CPx, NxCAx, _mm_fmadd_ps(CPy, NxCAy, _mm_mul_ps(CPz, NxCAz)));
  __m128 CP_dot_CA = _mm_fmadd_ps(CPx, CAx, _mm_fmadd_ps(CPy, CAy, _mm_mul_ps(CPz, CAz)));
  __m128 neg_CP_dot_BC = _mm_xor_ps(signbit, _mm_fmadd_ps(CPx, BCx, _mm_fmadd_ps(CPy, BCy, _mm_mul_ps(CPz, BCz))));

  __m128 out_ab = _mm_cmplt_ps(AP_dot_NxAB, _mm_setzero_ps());
  __m128 out_bc = _mm_andnot_ps(out_ab, _mm_cmplt_ps(BP_dot_NxBC, _mm_setzero_ps()));
  __m128 out_ab_bc = _mm_or_ps(out_ab, out_bc);
  __m128 out_ca = _mm_andnot_ps(out_ab_bc, _mm_cmplt_ps(CP_dot_NxCA, _mm_setzero_ps()));
  __m128 out_ab_bc_ca = _mm_or_ps(_mm_or_ps(out_ab, out_bc), out_ca);
  // out_none = !out_ab_bc_ca

  __m128 region = _mm_blendv_ps(_mm_blendv_ps(_mm_and_ps(_mm_set_ps1(1), out_ab),
                                              _mm_set_ps1(2), out_bc),
                                _mm_set_ps1(3), out_ca);

  __m128 w_a = _mm_blendv_ps(_mm_blendv_ps(_mm_andnot_ps(out_ab_bc_ca, BP_dot_NxBC),
                                           CP_dot_CA, out_ca),
                             neg_BP_dot_AB, out_ab);

  __m128 w_b = _mm_blendv_ps(_mm_blendv_ps(_mm_andnot_ps(out_ab_bc_ca, CP_dot_NxCA),
                                           neg_CP_dot_BC, out_bc),
                             AP_dot_AB, out_ab);

  __m128 w_c = _mm_blendv_ps(_mm_blendv_ps(_mm_andnot_ps(out_ab_bc_ca, AP_dot_NxAB),
                                           neg_AP_dot_CA, out_ca),
                             BP_dot_BC, out_bc);

  w_a = _mm_max_ps(w_a, _mm_setzero_ps());
  w_b = _mm_max_ps(w_b, _mm_setzero_ps());
  w_c = _mm_max_ps(w_c, _mm_setzero_ps());

  __m128 s = _mm_rcp_ps(_mm_add_ps(w_a, _mm_add_ps(w_b, w_c)));

  __m128 degenerate = _mm_cmpeq_ps(s, _mm_set1_ps(std::numeric_limits<float>::infinity()));
  if (_mm_movemask_ps(degenerate)) {  // probably quite unlikely, so it might beneficial to skip here. do profile.
    s = _mm_blendv_ps(s, _mm_set_ps1(1.f), degenerate);
    w_a = _mm_blendv_ps(s, _mm_set_ps1(1.f/3.f), degenerate);
    w_b = _mm_blendv_ps(s, _mm_set_ps1(1.f / 3.f), degenerate);
    w_c = _mm_blendv_ps(s, _mm_set_ps1(1.f / 3.f), degenerate);
  }

  __m128 Qx = _mm_mul_ps(s, _mm_fmadd_ps(w_a, Ax, _mm_fmadd_ps(w_b, Bx, _mm_mul_ps(w_c, Cx))));
  __m128 Qy = _mm_mul_ps(s, _mm_fmadd_ps(w_a, Ay, _mm_fmadd_ps(w_b, By, _mm_mul_ps(w_c, Cy))));
  __m128 Qz = _mm_mul_ps(s, _mm_fmadd_ps(w_a, Az, _mm_fmadd_ps(w_b, Bz, _mm_mul_ps(w_c, Cz))));

  Q_[0] = Qx.m128_f32[0];
  Q_[1] = Qy.m128_f32[0];
  Q_[2] = Qz.m128_f32[0];

  return region.m128_f32[0];
}


int main()
{
  const float V[][3] = 
  {
    {  0.2f,  0.15f, 0.0f },
    {  0.9f, -0.15f, 0.0f },
    { -0.9f, -0.15f, 0.0f },
  };

  const uint32_t ix[3] = { 0, 1, 2 };


  const float r = 2.f;
  const unsigned n = 40;
  const unsigned m = 20;


  std::ofstream mat("dump.mtl");
  for (unsigned i = 0; i < 8; i++) {
    mat << "newmtl m" << i << "\n";
    mat << "Kd "
      << (i & 1 ? 1 : 0) << ' '
      << (i & 2 ? 1 : 0) << ' '
      << (i & 4 ? 1 : 0) << '\n';
  }

  std::ofstream out("dump.obj");
  out << "mtllib dump.mtl\n";

  const auto & A = V[ix[0]];
  const auto & B = V[ix[1]];
  const auto & C = V[ix[2]];

  for (unsigned j = 0; j < m; j++) {
    for (unsigned i = 0; i < n; i++) {
      auto theta = (float(M_PI)*j + 0.5f) / m;
      auto phi = (float(2.0*M_PI)*i + 0.5f) / n;

      const float P[3] = {
        r * std::sin(theta)*std::cos(phi),
        r * std::sin(theta)*std::sin(phi),
        r * std::cos(theta)
      };
      float Q[3];
      auto region = nearestPointOnTriangle2(Q, A, B, C, P);
      //auto region = nearestPointOnTriangleSSE(Q, A, B, C, P);
      out << "usemtl m" << ((region + 1) & 7) << "\n";

      out << "v " << Q[0] << ' ' << Q[1] << ' ' << Q[2] << '\n';
      out << "v " << P[0] << ' ' << P[1] << ' ' << P[2] << '\n';
      out << "l -1 -2\n";

    }
  }

  // output triangle
  out << "usemtl m0\n";
  out << "v " << A[0] << ' ' << A[1] << ' ' << A[2] << '\n';
  out << "v " << B[0] << ' ' << B[1] << ' ' << B[2] << '\n';
  out << "v " << C[0] << ' ' << C[1] << ' ' << C[2] << '\n';
  out << "f -1 -2 -3\n";

  return 0;
}
