#ifndef _ALI_DCA_FITTER_
#define _ALI_DCA_FITTER_

#include <TMath.h>
#include <Rtypes.h>
#include "AliExternalTrackParam.h"

class DCAFitter {

public:
  
  // >>> Auxiliary structs used by DCA finder
  
  //----------------------------------------------------
  ///< Inverse cov matrix of the point defined by the track
  struct TrackPointCovI {
    double sxx, syy, syz, szz;
    
    TrackPointCovI(const AliExternalTrackParam& trc) { set(trc); }
    TrackPointCovI() {}
    void set(const AliExternalTrackParam& trc)
    {
      // we assign Y error to X for DCA calculation of 2 points
      // (otherwise for quazi-collinear tracks the X will not be constrained)
      double cyy = trc.GetSigmaY2(), czz = trc.GetSigmaZ2(), cyz = trc.GetSigmaZY(), cxx = cyy;
      double detYZ = cyy*czz - cyz*cyz;
      if (detYZ>0.) {
	sxx = 1./cxx;
	syy = czz/detYZ;
	syz = -cyz/detYZ;
	szz = cyy/detYZ;
      }
      else {
	syy = 0.0; // failure
      }
    }
    
    ClassDefNV(TrackPointCovI,1);
  };
  
  //----------------------------------------------------
  ///< particular point on track trajectory (in track proper alpha-frame)
  struct TrackPoint {
    double x, y, z;
    TrackPoint(AliExternalTrackParam& trc) { set(trc); }  
    TrackPoint(double px=0, double py=0, double pz=0) : x(px), y(py), z(pz) {}
    void set(const AliExternalTrackParam& trc) {
      x = trc.GetX();
      y = trc.GetY();
      z = trc.GetZ();
    }
    ClassDefNV(TrackPoint,1);
  };
  
  //----------------------------------------------------
  ///< Derivative (up to 2) of the TrackParam position over its running param X 
  struct TrackPointDeriv2 {
    double dydx, dzdx, d2ydx2, d2zdx2;
    
    TrackPointDeriv2() {}
    TrackPointDeriv2(const AliExternalTrackParam& trc, double bz) { set(trc,bz); }
    void set(const AliExternalTrackParam& trc, double bz) {
      double snp = trc.GetSnp(), csp = TMath::Sqrt((1.-snp)*(1.+snp)), cspI = 1./csp, crv2c = trc.GetC(bz)*cspI;
      dydx = snp*cspI;           // = snp/csp
      dzdx = trc.GetTgl()*cspI; // = tgl/csp
      d2ydx2 = crv2c*cspI*cspI;  // = crv/csp^3
      d2zdx2 = crv2c*dzdx*dydx;  // = crv*tgl*snp/csp^3
    }
    
    ClassDefNV(TrackPointDeriv2,1);
  };
  
  //----------------------------------------------------
  //< precalculated track radius, center, alpha sin,cos and their combinations
  struct TrackAuxPar {
    double c,s; // cos ans sin of track alpha
    double cc, cs, ss; // products
    double r, xCen, yCen; // helix radius and center in lab
    
    TrackAuxPar() {}
    TrackAuxPar(const AliExternalTrackParam& trc, double bz) { set(trc,bz); }

    void set(const AliExternalTrackParam& trc, double bz) {
      c = TMath::Cos(trc.GetAlpha());
      s = TMath::Sin(trc.GetAlpha());
      setRCen(trc, bz);
      cc = c*c;
      ss = s*s;
      cs = c*s;
    }
    
    void setRCen(const AliExternalTrackParam& tr, double bz);
    
    void glo2loc(double vX,double vY, double& vXL, double& vYL) const {
      // rotate XY in global frame to the frame of track with angle A
      vXL = vX*c + vY*s;
      vYL = -vX*s + vY*c;
    }
    
    void loc2glo(double vXL,double vYL, double& vX, double& vY) const {
      // rotate XY in local alpha frame to global frame
      vX = vXL*c - vYL*s;
      vY = vXL*s + vYL*c;
    }

    ClassDefNV(TrackAuxPar,1);
  };
  
  //----------------------------------------------------
  //< crossing coordinates of 2 circles
  struct CrossInfo {
    double xDCA[2];
    double yDCA[2];
    int nDCA;
    
    CrossInfo() {}
    CrossInfo(const TrackAuxPar& trc0, const TrackAuxPar& trc1) { set(trc0, trc1); }
    void set(const TrackAuxPar& trc0, const TrackAuxPar& trc1);
    
    ClassDefNV(CrossInfo,1);
  };

  struct Derivatives {
    double dChidx0, dChidx1; // 1st derivatives of chi2 vs tracks local parameters X
    double dChidx0dx0, dChidx1dx1, dChidx0dx1; // 2nd derivatives of chi2 vs tracks local parameters X
  };

  //----------------------------------------------------
  //< coefficients of the track-point contribution to the PCA (Vx,Vy,Vz) to 2 points in lab frame represented via local points coordinates as
  //< Vx = mXX0*x0+mXY0*y0+mXZ0*z0 + mXX1*x1+mXY1*y1+mXZ1*z1
  //< Vy = mYX0*x0+mYY0*y0+mYZ0*z0 + mYX1*x1+mYY1*y1+mYZ1*z1
  //< Vz = mZX0*x0+mZY0*y0+mZZ0*z0 + mZX1*x1+mZY1*y1+mZZ1*z1
  //< where {x0,y0,z0} and {x1,y1,z1} are track positions in their local frames
  struct TrackCoefVtx {
    double mXX,mXY,mXZ,mYX,mYY,mYZ,mZX,mZY,mZZ;
    TrackCoefVtx() {}
    ClassDefNV(TrackCoefVtx,1);
  };
  
  // <<< Auxiliary structs used by DCA finder

  //===============================================================================
  
  DCAFitter() {};

  double getMaxIter() const { return mMaxIter; }
  double getMaxR()    const { return TMath::Sqrt(mMaxR2); }
  double getMaxChi2() const { return mMaxChi2;}
  double getMinParamChange() const { return mMinParamChange; }
  double getBz() const { return mBz; }
  bool   getUseAbsDCA() const {return mUseAbsDCA;}
  
  void setMaxIter(int n=20) { mMaxIter = n>2 ? n : 2; }
  void setMaxR(double r=200.) { mMaxR2 = r*r; }
  void setMaxChi2(double chi2=999.) { mMaxChi2 = chi2; }
  void setBz(double bz) { mBz = bz; }
  void setMinParamChange(double x=1e-3) { mMinParamChange = x>1e-4 ? x : 1.e-4; }
  void setMinRelChi2Change(double r=0.9) { mMinRelChi2Change = r>0.1 ? r : 999.;}
  void setUseAbsDCA(bool v) {mUseAbsDCA = v;}
    
  DCAFitter(double bz, double minRelChiChange=0.9, double minXChange=1e-3, double maxChi=999, int n=20, double maxR = 200.) {
    setMaxIter(n);
    setMaxR(maxR);
    setMaxChi2(maxChi);
    setMinParamChange(minXChange);
    setMinRelChi2Change(minRelChiChange);
    setBz(bz);
    setUseAbsDCA(false); // by default use weighted DCA definition (much slower)
  }

  ///< number of validated V0 candidates (at most 2 are possible)
  int getNCandidates() const { return mNCandidates; }

  ///< return PCA candidate (no check for its validity)
  TrackPoint& getPCACandidate(int cand) { return mPCA[cand]; }

  ///< return Chi2 at PCA candidate (no check for its validity)
  double getChi2AtPCACandidate(int cand) { return mChi2[cand]; }

  ///< 1st track params propagated to V0 candidate (no check for the candidate validity)
  AliExternalTrackParam& getTrack0(int cand) { return mCandTr0[cand]; }

  ///< 2nd track params propagated to V0 candidate (no check for the candidate validity)
  AliExternalTrackParam& getTrack1(int cand) { return mCandTr1[cand]; }  

  
  ///< calculate parameters tracks at PCA
  int process(const AliExternalTrackParam& trc0, const AliExternalTrackParam& trc1);

  ///< calculate parameters tracks at PCA, using precalculated aux info // = TrackAuxPar(track) //
  int process(const AliExternalTrackParam& trc0, const TrackAuxPar& trc0Aux,
	      const AliExternalTrackParam& trc1, const TrackAuxPar& trc1Aux);

  ///< minimizer for abs distance definition of DCA
  bool processCandidateDCA(const TrackAuxPar& trc0Aux, const TrackAuxPar& trc1Aux);

  ///< minimizer for weighted distance definition of DCA (chi2)
  bool processCandidateChi2(const TrackAuxPar& trc0Aux, const TrackAuxPar& trc1Aux);


 protected:
  void calcPCACoefs(const TrackAuxPar& trc0Aux, const TrackPointCovI& trcEI0,
		    const TrackAuxPar& trc1Aux, const TrackPointCovI& trcEI1,
		    TrackCoefVtx& trCFVT0, TrackCoefVtx& trCFVT1) const;

  ///< PCA with weighted DCA definition
  void calcPCA(const AliExternalTrackParam& trc0, const TrackCoefVtx& trCFVT0,
	       const AliExternalTrackParam& trc1, const TrackCoefVtx& trCFVT1,
	       double &Vx, double &Vy, double &Vz) const;

  ///< PCA with weighted DCA definition
  void calcPCA(const TrackPoint& tPnt0, const TrackCoefVtx& trCFVT0,
	       const TrackPoint& tPnt1, const TrackCoefVtx& trCFVT1,
	       double &Vx, double &Vy, double &Vz) const;

  ///< PCA with abs DCA definition
  void calcPCA(const AliExternalTrackParam& trc0, const TrackAuxPar& trc0Aux, 
	       const AliExternalTrackParam& trc1, const TrackAuxPar& trc1Aux, 
	       double &Vx, double &Vy, double &Vz) const;
  
  ///< PCA with abs DCA definition
  void calcPCA(const TrackPoint& tPnt0, const TrackAuxPar& trc0Aux, 
	       const TrackPoint& tPnt1, const TrackAuxPar& trc1Aux, 
	       double &Vx, double &Vy, double &Vz) const;

  ///< chi2 (weighted distance)
  double calcChi2(double Vx, double Vy, double Vz,
		  const TrackPoint& tPnt0, const TrackAuxPar& trc0Aux, const TrackPointCovI& trcEI0,
		  const TrackPoint& tPnt1, const TrackAuxPar& trc1Aux, const TrackPointCovI& trcEI1) const;

  ///< DCA (abs distance)
  double calcDCA(const TrackPoint& tPnt0, const TrackAuxPar& trc0Aux, const TrackPoint& tPnt1, const TrackAuxPar& trc1Aux) const;
  double calcDCA(double Vx, double Vy, double Vz, const TrackPoint& tPnt0, const TrackAuxPar& trc0Aux, const TrackPoint& tPnt1, const TrackAuxPar& trc1Aux) const;
  
  TrackPoint calcResid(const TrackPoint& pnt, const TrackAuxPar& alpCS, const TrackPoint& vtx) const
  {
    double vlX,vlY; // Vertex XY in track local frame
    alpCS.glo2loc(vtx.x,vtx.y,vlX,vlY);
    return TrackPoint(pnt.x - vlX, pnt.y - vlY, pnt.z - vtx.z);
  }

  void chi2Deriv(const TrackPoint& tPnt0, const TrackPointDeriv2& tDer0, const TrackAuxPar& trc0Aux, const TrackPointCovI& trcEI0, const TrackCoefVtx& trCFVT0,
		 const TrackPoint& tPnt1, const TrackPointDeriv2& tDer1, const TrackAuxPar& trc1Aux, const TrackPointCovI& trcEI1, const TrackCoefVtx& trCFVT1,
		 Derivatives& deriv) const;

  void DCADeriv(const TrackPoint& tPnt0, const TrackPointDeriv2& tDer0, const TrackAuxPar& trc0Aux, 
		const TrackPoint& tPnt1, const TrackPointDeriv2& tDer1, const TrackAuxPar& trc1Aux, 
		Derivatives& deriv) const;

  bool closerToAlternative(float x, float y) const;
  
 private:
  bool   mUseAbsDCA; // ignore track errors (minimize abs DCA to vertex)
  double mMaxIter; // max iterations
  double mMaxR2;   // max radius to consider 
  double mMaxChi2; // max chi2 to accept
  double mMinRelChi2Change; // stop iterations if relative chi2 change is less than requested
  double mMinParamChange; // stop iterations when both X params change by less than this value

  double mBz;      // mag field for simple propagation

  CrossInfo mCrossings; // analystical XY crossings (max 2) of the seeds
  int mCrossITCur;   // XY crossing being tested
  int mCrossIDAlt;   // XY crossing alternative to the one being tested. Abandon fit if it converges to it
  
  int mNCandidates; // number of consdered candidates
  AliExternalTrackParam mCandTr0[2], mCandTr1[2]; // Tracks at PCA: max 2 candidates possible
  TrackPoint mPCA[2]; // PCA for 2 possible cases
  double mChi2[2]; // Chi2 at PCA candidate
  TrackPointCovI mTrcEI0[2], mTrcEI1[2]; // errors for each track candidate
  TrackCoefVtx mTrCFVT0[2], mTrCFVT1[2]; // coefficients of PCA vs track points for each track
  
  ClassDefNV(DCAFitter,1);
};


#endif
