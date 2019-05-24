#include "AliExternalTrackParam.h"

const double kTinyD = 1e-12;


///< Inverse cov matrix of the point defined by the track
struct TrackPointCovI {
  double sxx, syy, syz, szz;

  void set(const AliExternalTrackParam& trc)
  {
    // we assign Y error to X for DCA calculation of 2 points
    // (otherwise for quazi-collinear tracks the X will not be constrained)
    double cyy = trc.GetSigmaY2(), czz = trc.GetSigmaZ2(), cyz = trc.GetSigmaZY(), cxx = cyy;
    double detYZ = cyy*czz - cyz*cyz;
    if (detYZ>kTinyD) {
      sxx = 1./cxx;
      syy = czz/detYZ;
      syz = -cyz/detYZ;
      szz = cyy/detYZ;
    }
    else {
      syy = 0.0; // failure
    }
  }
  TrackPointCovI(const AliExternalTrackParam& trc)
  {
    set(trc);
  }
  TrackPointCovI() {}

  ClassDefNV(TrackPointCovI,1);
};

///< particular point on track trajectory (in track proper alpha-frame)
struct TrackPoint {
  double x, y, z;

  TrackPoint(double px=0, double py=0, double pz=0) : x(px), y(py), z(pz) {}
  TrackPoint(AliExternalTrackParam& trc) { set(trc); }
  void set(const AliExternalTrackParam& trc) {
    x = trc.GetX();
    y = trc.GetY();
    z = trc.GetZ();
  }
  ClassDefNV(TrackPoint,1);
};

///< Derivative (up to 2) of the TrackParam position over its running param X 
struct TrackPointDeriv2 {

  double dydx, dzdx, d2ydx2, d2zdx2;
  void setDeriv(const AliExternalTrackParam& trc, double bz)
  {
    double snp = trc.GetSnp(), csp = TMath::Sqrt((1.-snp)*(1.+snp)), cspI = 1./csp, crv2c = trc.GetC(bz)*cspI;
    dydx = snp*cspI;           // = snp/csp
    dzdx = trc.GetTgl()*cspI; // = tgl/csp
    d2ydx2 = crv2c*cspI*cspI;  // = crv/csp^3
    d2zdx2 = crv2c*dzdx*dydx;  // = crv*tgl*snp/csp^3
  }

  TrackPointDeriv2() {}
  TrackPointDeriv2(const AliExternalTrackParam& trc, double bz) {
    setDeriv(trc,bz);
  }
  
  ClassDefNV(TrackPointDeriv2,1);
};

///< Derivative (up to 3) of the TrackParam position over its running param X 
struct TrackPointDeriv3 {

  double dydx, dzdx, d2ydx2, d2zdx2, d3ydx3, d3zdx3; // make sure we need this
  void setDeriv(const AliExternalTrackParam& trc, double bz)
  {
    double snp = trc.GetSnp(), csp = TMath::Sqrt((1.-snp)*(1.+snp)), cspI = 1./csp, crv2c = trc.GetC(bz)*cspI;
    dydx = snp*cspI;           // = snp/csp
    dzdx = trc.GetTgl()*cspI; // = tgl/csp
    d2ydx2 = crv2c*cspI*cspI;  // = crv/csp^3
    d2zdx2 = crv2c*dzdx*dydx;  // = crv*tgl*snp/csp^3
    d3ydx3 = 3.*crv2c*d2ydx2*dydx; // 3*crv^2*snp/csp^5
    d3zdx3 = 3.*crv2c*d2zdx2*dydx; // 3*crv^2*tgl*snp^2/crv^5
  }

  TrackPointDeriv3() {}
  TrackPointDeriv3(const AliExternalTrackParam& trc, double bz) {
    setDeriv(trc,bz);
  }
  
  ClassDefNV(TrackPointDeriv3,1);
};

//< precalculate track alpha sin,cos and combinations
struct TrackAlpSinCos {
  double c,s; // cos ans sin of track alpha
  double cc, cs, ss; // products
  void set(double alp) {
    c = TMath::Cos(alp);
    s = TMath::Sin(alp);
    cc = c*c;
    ss = s*s;
    cs = c*s;
  }

  TrackAlpSinCos() {}
  TrackAlpSinCos(double alp) {
    set(alp);
  }

  void glo2local(double vX,double vY, double& vXL, double& vYL) const {
    // rotate XY in global frame to the frame of track with angle A
    vXL = vX*c + vY*s;
    vYL = -vX*s + vY*c;
  }
  ClassDefNV(TrackAlpSinCos,1);
};


//< coefficients of the track-point contribution to the PCA (Vx,Vy,Vz) to 2 points in lab frame represented via local points coordinates as
//< Vx = mXX0*x0+mXY0*y0+mXZ0*z0 + mXX1*x1+mXY1*y1+mXZ1*z1
//< Vy = mYX0*x0+mYY0*y0+mYZ0*z0 + mYX1*x1+mYY1*y1+mYZ1*z1
//< Vz = mZX0*x0+mZY0*y0+mZZ0*z0 + mZX1*x1+mZY1*y1+mZZ1*z1
//< where {x0,y0,z0} and {x1,y1,z1} are track positions in their local frames
struct TrackCoefVtx {
  double mXX,mXY,mXZ,mYX,mYY,mYZ,mZX,mZY,mZZ;
  ClassDefNV(TrackCoefVtx,1);
};


struct DCAFitter {

  TrackPoint       mP0, mP1;     // track poitns where current PCA is evaluated
  TrackAlpSinCos   mCS0, mCS1;   // precalculated tracks alpha frame sin,cos ...
  TrackPointDeriv2 mD0, mD1;     // derivatives of the tracks vs X params
  TrackPointCovI   mEI0, mEI1;   // inverted error matrices of the track points
  TrackCoefVtx     mVT0, mVT1;   // coefficients of the track-point contribution to the PCA 

  DCAFitter() {};
  void calcPCACoefs();
  void calcPCA(double &Vx, double &Vy, double &Vz);
  void calcPCA(AliExternalTrackParam& trc0, AliExternalTrackParam& trc1, double &Vx, double &Vy, double &Vz);

  TrackPoint calcResid(const TrackPoint& pnt, const TrackAlpSinCos& alpCS, const TrackPoint& vtx) const {
    double vlX,vlY; // Vertex XY in track local frame
    alpCS.glo2local(vtx.x,vtx.y,vlX,vlY);
    return TrackPoint(pnt.x - vlX, pnt.y - vlY, pnt.z - vtx.z);
  }

  void chi2Deriv(AliExternalTrackParam& trc0, AliExternalTrackParam& trc1, double bz);
  
};


void DCAFitter::calcPCA(AliExternalTrackParam& trc0, AliExternalTrackParam& trc1, double &Vx, double &Vy, double &Vz)
{
  // calculate coefficients of the PCA (Vx,Vy,Vz) to 2 points in lab frame represented via local points coordinates
  Vx = mVT0.mXX*trc0.GetX() + mVT0.mXY*trc0.GetY() + mVT0.mXZ*trc0.GetX()
    +  mVT1.mXX*trc1.GetX() + mVT1.mXY*trc1.GetY() + mVT1.mXZ*trc1.GetX();
  Vx = mVT0.mYX*trc0.GetX() + mVT0.mYY*trc0.GetY() + mVT0.mYZ*trc0.GetX()
    +  mVT1.mYX*trc1.GetX() + mVT1.mYY*trc1.GetY() + mVT1.mYZ*trc1.GetX();
  Vx = mVT0.mZX*trc0.GetX() + mVT0.mZY*trc0.GetY() + mVT0.mZZ*trc0.GetX()
    +  mVT1.mZX*trc1.GetX() + mVT1.mZY*trc1.GetY() + mVT1.mZZ*trc1.GetX();
}

void DCAFitter::calcPCA(double &Vx, double &Vy, double &Vz)
{
  // calculate coefficients of the PCA (Vx,Vy,Vz) to 2 points in lab frame represented via local points coordinates
  Vx = mVT0.mXX*mP0.x + mVT0.mXY*mP0.y + mVT0.mXZ*mP0.z
    +  mVT1.mXX*mP1.x + mVT1.mXY*mP1.y + mVT1.mXZ*mP1.z;
  Vx = mVT0.mYX*mP0.x + mVT0.mYY*mP0.y + mVT0.mYZ*mP0.z
    +  mVT1.mYX*mP1.x + mVT1.mYY*mP1.y + mVT1.mYZ*mP1.z;
  Vx = mVT0.mZX*mP0.x + mVT0.mZY*mP0.y + mVT0.mZZ*mP0.z
    +  mVT1.mZX*mP1.x + mVT1.mZY*mP1.y + mVT1.mZZ*mP1.z;
}

void DCAFitter::chi2Deriv(AliExternalTrackParam& trc0, AliExternalTrackParam& trc1, double bz)
{

  //mD0.setDeriv2(trc0, bz); // 1st and 2nd derivatives of 1st track over its X param
  //mD1.setDeriv2(trc1, bz); // 1st and 2nd derivatives of 2nd track over its X param

  // mP0.set(trc0);
  // mP1.set(trc1);
  
  TrackPoint vtx;
  calcPCA(vtx.x,vtx.y,vtx.z); // calculate PCA for current track-points positions
  // calculate residuals
  auto res0 = calcResid(mP0,mCS0,vtx);
  auto res1 = calcResid(mP1,mCS1,vtx);
  
  // res0.x,dy0,dz0 = x0 - Xl0, y0 - Yl0, ... , whith x0,y0,z0: track coords in its frame, and Xl0, Yl0, Zl0 : vertex rotated to same frame
  
  // aux params to minimize multiplications  
  double dx0s = res0.x*mEI0.sxx;
  double dyz0sY = res0.y*mEI0.syy + res0.z*mEI0.syz;
  double dyz0sZ = res0.y*mEI0.syz + res0.z*mEI0.szz;
  double dx1s = res1.x*mEI1.sxx;
  double dyz1sY = res1.y*mEI1.syy + res1.z*mEI1.syz;
  double dyz1sZ = res1.y*mEI1.syz + res1.z*mEI1.szz;

  double xx0DtXYXZx0 = mVT0.mXX + mVT0.mXY*mD0.dydx + mVT0.mXZ*mD0.dzdx;
  double yz0DtYYYZx0 = mVT0.mYX + mVT0.mYY*mD0.dydx + mVT0.mYZ*mD0.dzdx;
  double xx1DtXYXZx1 = mVT1.mXX + mVT1.mXY*mD1.dydx + mVT1.mXZ*mD1.dzdx;
  double yz1DtYYYZx1 = mVT1.mYX + mVT1.mYY*mD1.dydx + mVT1.mYZ*mD1.dzdx;
  
  double DtXYXZx02 = mVT0.mXY*mD0.d2ydx2 + mVT0.mXZ*mD0.d2zdx2;
  double DtYYYZx02 = mVT0.mYY*mD0.d2ydx2 + mVT0.mYZ*mD0.d2zdx2;
  double DtXYXZx12 = mVT1.mXY*mD1.d2ydx2 + mVT1.mXZ*mD1.d2zdx2;
  double DtYYYZx12 = mVT1.mYY*mD1.d2ydx2 + mVT1.mYZ*mD1.d2zdx2;

  double FDdx0Dx0 = 1. - mCS0.c*xx0DtXYXZx0 - mCS0.s*yz0DtYYYZx0;
  double FDdx1Dx1 = 1. - mCS1.c*xx1DtXYXZx1 - mCS1.s*yz1DtYYYZx1;

  double FDdy0Dx0 = mD0.dydx + mCS0.s*xx0DtXYXZx0 - mCS0.c*yz0DtYYYZx0;
  double FDdy1Dx1 = mD1.dydx + mCS1.s*xx1DtXYXZx1 - mCS1.c*yz1DtYYYZx1;

  double FDdz0Dx0 = -mVT0.mZX - mVT0.mZY*mD0.dydx + mD0.dzdx - mVT0.mZZ*mD0.dzdx;
  double FDdz1Dx1 = -mVT1.mZX - mVT1.mZY*mD1.dydx + mD1.dzdx - mVT1.mZZ*mD1.dzdx;

  double FDdx0Dx0x0 = -(mCS0.c*DtXYXZx02 + mCS0.s*DtYYYZx02);
  double FDdy0Dx0x0 =  (mCS0.s*DtXYXZx02 - mCS0.c*DtYYYZx02) + mD0.d2ydx2;
  double FDdz0Dx0x0 = -mVT0.mZY*mD0.d2ydx2 + mD0.d2zdx2*(1.-mVT0.mZZ);

  double FDdx1Dx1x1 = -(mCS1.c*DtXYXZx12 + mCS1.s*DtYYYZx12);
  double FDdy1Dx1x1 =  (mCS1.s*DtXYXZx12 - mCS1.c*DtYYYZx12) + mD1.d2ydx2;
  double FDdz1Dx1x1 = -mVT1.mZY*mD1.d2ydx2 + mD1.d2zdx2*(1.-mVT1.mZZ);

  // aux params to minimize multiplications
  double FD0YYYZ =  FDdy0Dx0*mEI0.syy + FDdz0Dx0*mEI0.syz;
  double FD0YZZZ =  FDdy0Dx0*mEI0.syz + FDdz0Dx0*mEI0.szz;
  double FD1YYYZ =  FDdy1Dx1*mEI1.syy + FDdz1Dx1*mEI1.syz;
  double FD1YZZZ =  FDdy1Dx1*mEI1.syz + FDdz1Dx1*mEI1.szz;
  
  // 1st derivatives over track params x
  double dChidx0 = dx0s*FDdx0Dx0 + res0.y*FD0YYYZ + res0.z*FD0YZZZ;
  double dChidx1 = dx1s*FDdx1Dx1 + res1.y*FD1YYYZ + res1.z*FD0YZZZ;

  // 2nd derivative over track params x
  double dChidx0dx0 =
    FDdx0Dx0x0*dx0s +
    FDdy0Dx0x0*dyz0sY +
    FDdz0Dx0x0*dyz0sZ +    
    FDdx0Dx0*FDdx0Dx0*mEI0.sxx +
    FDdy0Dx0*FD0YYYZ + 
    FDdz0Dx0*FD0YZZZ;

  double dChidx1dx1 =
    FDdx1Dx1x1*dx1s +
    FDdy1Dx1x1*dyz1sY +
    FDdz1Dx1x1*dyz1sZ +    
    FDdx1Dx1*FDdx1Dx1*mEI1.sxx +
    FDdy1Dx1*FD1YYYZ + 
    FDdz1Dx1*FD1YZZZ;
}


void DCAFitter::calcPCACoefs()
{
  // calculate coefficients of the PCA (Vx,Vy,Vz) to 2 points in lab frame represented via local points coordinates as
  // Vx = mXX0*x0+mXY0*y0+mXZ0*z0 + mXX1*x1+mXY1*y1+mXZ1*z1
  // Vy = mYX0*x0+mYY0*y0+mYZ0*z0 + mYX1*x1+mYY1*y1+mYZ1*z1
  // Vz = mZX0*x0+mZY0*y0+mZZ0*z0 + mZX1*x1+mZY1*y1+mZZ1*z1
  // where {x0,y0,z0} and {x1,y1,z1} are track positions in their local frames
  //
  // we find the PCA of 2 tracks poins weighted by their errors, i.e. minimizing
  // chi2 = ....
  // these are the coefficients of dChi2/d{Vx,Vy,vZ} = 0
  double axx = mCS0.cc*mEI0.sxx + mCS1.cc*mEI1.sxx + mCS0.ss*mEI0.syy + mCS1.ss*mEI1.syy;
  double axy = mCS0.cs*(mEI0.sxx - mEI0.syy) + mCS1.cs*(mEI1.sxx - mEI1.syy);
  double axz = -(mCS0.s*mEI0.syz + mCS1.s*mEI1.syz);
  double ayy = mCS0.ss*mEI0.sxx + mCS1.ss*mEI1.sxx + mCS0.cc*mEI0.syy + mCS1.cc*mEI1.syy; // = (mEI0.sxx + mEI1.sxx + mEI0.syy + mEI1.syy) - axx
  double ayz = mCS0.c*mEI0.syz + mCS1.c*mEI1.syz;
  double azz = mEI0.szz + mEI1.szz;
  //
  // define some aux variables
  double axxyy = axx*ayy, axxzz = axx*azz, axxyz = axx*ayz, 
    axyxy = axy*axy, axyxz = axy*axz, axyyz = axy*ayz, axyzz = axy*azz, 
    axzxz = axz*axz, axzyy = axz*ayy, axzyz = axz*ayz,
    ayzyz = ayz*ayz,  ayyzz = ayy*azz;
  double dAxxyyAxyxy = axxyy - axyxy, dAxyyzAxzyy = axyyz - axzyy, dAxyxzAxxyz = axyxz - axxyz;
  double dAxzyzAxyzz = axzyz - axyzz, dAyyzzAyzyz = ayzyz - ayyzz, dAxxzzAxzxz = axxzz - axzxz;
  double det = -dAxyyzAxzyy*axz + dAxyxzAxxyz*ayz + dAxxyyAxyxy*azz, detI = 1./det;
  
  double dfxPCS0 = dAyyzzAyzyz*mCS0.c + dAxzyzAxyzz*mCS0.s, dfxQCS0 = dAxzyzAxyzz*mCS0.c + dAyyzzAyzyz*mCS0.s;
  double dfyPCS0 = dAxzyzAxyzz*mCS0.c + dAxxzzAxzxz*mCS0.s, dfyQCS0 = dAxxzzAxzxz*mCS0.c - dAxzyzAxyzz*mCS0.s;
  double dfzPCS0 = dAxyyzAxzyy*mCS0.c + dAxyxzAxxyz*mCS0.s, dfzQCS0 = dAxyxzAxxyz*mCS0.c - dAxyyzAxzyy*mCS0.s;

  double dfxPCS1 = dAyyzzAyzyz*mCS1.c + dAxzyzAxyzz*mCS1.s, dfxQCS1 = dAxzyzAxyzz*mCS1.c + dAyyzzAyzyz*mCS1.s;
  double dfyPCS1 = dAxzyzAxyzz*mCS1.c + dAxxzzAxzxz*mCS1.s, dfyQCS1 = dAxxzzAxzxz*mCS1.c - dAxzyzAxyzz*mCS1.s;
  double dfzPCS1 = dAxyyzAxzyy*mCS1.c + dAxyxzAxxyz*mCS1.s, dfzQCS1 = dAxyxzAxxyz*mCS1.c - dAxyyzAxzyy*mCS1.s;
  //
  mVT0.mXX = detI*(dfxPCS0*mEI0.sxx);
  mVT0.mXY = detI*(dfxQCS0*mEI0.syy + dAxyyzAxzyy*mEI0.syz);
  mVT0.mXZ = detI*(dfxQCS0*mEI0.syz + dAxyyzAxzyy*mEI0.szz);
  
  mVT0.mYX = detI*(dfyPCS0*mEI0.sxx);
  mVT0.mYY = detI*(dfyQCS0*mEI0.syy + dAxyxzAxxyz*mEI0.syz);
  mVT0.mYZ = detI*(dfyQCS0*mEI0.syz + dAxyxzAxxyz*mEI0.szz);
  
  mVT0.mZX = detI*(dfzPCS0*mEI0.sxx);
  mVT0.mZY = detI*(dfzQCS0*mEI0.syy + dAxxyyAxyxy*mEI0.syz);
  mVT0.mZZ = detI*(dfzQCS0*mEI0.syz + dAxxyyAxyxy*mEI0.szz);
  
  mVT0.mXX = detI*(dfxPCS1*mEI1.sxx);
  mVT0.mXY = detI*(dfxQCS1*mEI1.syy + dAxyyzAxzyy*mEI1.syz);
  mVT0.mXZ = detI*(dfxQCS1*mEI1.syz + dAxyyzAxzyy*mEI1.szz);
  
  mVT0.mYX = detI*(dfyPCS1*mEI1.sxx);
  mVT0.mYY = detI*(dfyQCS1*mEI1.syy + dAxyxzAxxyz*mEI1.syz);
  mVT0.mYZ = detI*(dfyQCS1*mEI1.syz + dAxyxzAxxyz*mEI1.szz);
  
  mVT0.mZX = detI*(dfzPCS1*mEI1.sxx);
  mVT0.mZY = detI*(dfzQCS1*mEI1.syy + dAxxyyAxyxy*mEI1.syz);
  mVT0.mZZ = detI*(dfzQCS1*mEI1.syz + dAxxyyAxyxy*mEI1.szz);

}

