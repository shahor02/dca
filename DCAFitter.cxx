#include "DCAFitter.h"

const double kTinyD = 1e-12;
const double kDummyVal = 1e9;

//_____________________________________________________________________________________
void DCAFitter::CrossInfo::set(const TrackAuxPar& trc0, const TrackAuxPar& trc1)
{
  // calculate up to 2 crossings between 2 circles
  nDCA = 0;
  const auto& trcA = trc0.r>trc1.r ? trc0 : trc1;  // designate the largest circle as A
  const auto& trcB = trc0.r>trc1.r ? trc1 : trc0;
  double xDist = trcB.xCen - trcA.xCen, yDist = trcB.yCen - trcA.yCen;
  double dist2 = xDist*xDist + yDist*yDist, dist = TMath::Sqrt(dist2), rsum = trcA.r+trcB.r;
  if (TMath::Abs(dist)<kTinyD) return; // circles are concentric?
  double distI = 1./dist;
  if (dist > rsum) { // circles don't touch, chose a point in between
    // the parametric equation of lines connecting the centers is
    // x = x0 + t/dist * (x1-x0), y = y0 + t/dist * (y1-y0)
    nDCA = 1;
    xDCA[0] = 0.5*(trcA.xCen + xDist*distI*(dist+trcA.r-trcB.r));
    yDCA[0] = 0.5*(trcA.yCen + yDist*distI*(dist+trcA.r-trcB.r));
  }
  else if (dist+trcB.r < trcA.r) { // the small circle is nestled into large one w/o touching
    // select the point of closest approach of 2 circles
    nDCA = 1;
    xDCA[0] = 0.5*(trcB.xCen + trcA.xCen + xDist*distI*rsum);
    yDCA[0] = 0.5*(trcB.yCen + trcA.yCen + yDist*distI*rsum);
  }
  else { // 2 intersection points
    // to simplify calculations, we move to new frame x->x+Xc0, y->y+Yc0, so that
    // the 1st one is centered in origin
    if (TMath::Abs(xDist)<TMath::Abs(yDist)) {
      double a = (trcA.r*trcA.r-trcB.r*trcB.r+dist2)/(2.*yDist), b = -xDist/yDist, ab=a*b,bb=b*b;
      double det = ab*ab - (1.+bb)*(a*a-trcA.r*trcA.r);
      if (det>=0.) {
	if (det>0.) {
	  det = TMath::Sqrt(det);
	  xDCA[0] = (-ab + det)/(1.+b*b);
	  yDCA[0] = a + b*xDCA[0] + trcA.yCen;
	  xDCA[0] += trcA.xCen;
	  xDCA[1] = (-ab - det)/(1.+b*b);
	  yDCA[1] = a + b*xDCA[1] + trcA.yCen;
	  xDCA[1] += trcA.xCen;
	  nDCA = 2;
	}
	else {
	  xDCA[0] = -ab/(1.+bb);
	  yDCA[0] = a + b*xDCA[0] + trcA.yCen;
	  xDCA[0] += trcA.xCen;
	  nDCA = 1;
	}
      }
      else {
	printf("This should not happen\n");
	return;
      }
    }
    else {
      double a = (trcA.r*trcA.r-trcB.r*trcB.r+dist2)/(2.*xDist), b = -yDist/xDist, ab=a*b,bb=b*b;
      double det = ab*ab - (1.+bb)*(a*a-trcA.r*trcA.r);
      if (det>=0.) {
	if (det>0.) {
	  det = TMath::Sqrt(det);
	  yDCA[0] = (-ab + det)/(1.+bb);
	  xDCA[0] = a + b*yDCA[0] + trcA.xCen;
	  yDCA[0] += trcA.yCen;
	  yDCA[1] = (-ab - det)/(1.+bb);
	  xDCA[1] = a + b*yDCA[1] + trcA.xCen;
	  yDCA[1] += trcA.yCen;
	  nDCA = 2;
	}
	else {
	  yDCA[0] = -ab/(1.+bb);
	  xDCA[0] = a + b*yDCA[0] + trcA.xCen;
	  yDCA[0] += trcA.yCen;
	    nDCA = 1;
	}
      }
      else {
	printf("This should not happen\n");
	return;
      }
    }
  }  
}

//_____________________________________________________________________________________
void DCAFitter::TrackAuxPar::setRCen(const AliExternalTrackParam& tr, double bz)
{
  // set track radius and circle coordinates in global frame
  double crv = tr.GetC(bz);
  r = TMath::Abs(crv);
  r = 1./r;
  double sn = tr.GetSnp();
  double cs = TMath::Sqrt( (1.-sn)*(1.+sn));
  double x,y;
  if (crv>0) { // clockwise
    x = tr.GetX() - sn*r;
    y = tr.GetY() + cs*r;
  }
  else {
    x = tr.GetX() + sn*r;
    y = tr.GetY() - cs*r;    
  }
  loc2glo(x,y,xCen,yCen);
}

//___________________________________________________________________
int DCAFitter::process(const AliExternalTrackParam& trc0, const AliExternalTrackParam& trc1)
{
  // find dca of 2 tracks
  TrackAuxPar trc0Aux(trc0,mBz), trc1Aux(trc1,mBz);
  return process(trc0, trc0Aux, trc1, trc1Aux);
}

//___________________________________________________________________
int DCAFitter::process(const AliExternalTrackParam& trc0, const TrackAuxPar& trc0Aux,
		       const AliExternalTrackParam& trc1, const TrackAuxPar& trc1Aux)
{
  // find dca of 2 tracks with aux info preclaculated
  mNCandidates = 0;
  mCrossings.set(trc0Aux, trc1Aux); // find at most 2 candidates of 2 circles crossing
  if (!mCrossings.nDCA) {
    return 0; // no crossing
  }
  for (int ic=0;ic<mCrossings.nDCA;ic++) {
    // both XY crossings may eventually converge to the same point. To stop asap this redundant step
    // we abandon fit if it appears to be closer to alternative crossing
    mCrossITCur = ic;
    mCrossIDAlt = mCrossings.nDCA==2 ? 1-ic : -1;
    // check if radius is acceptable
    if (mCrossings.xDCA[ic]*mCrossings.xDCA[ic]+mCrossings.yDCA[ic]*mCrossings.yDCA[ic] > mMaxR2) {
      continue;
    }
    // find dca starting from proximity of transverse point xv,yv
    double xl,yl;
    mCandTr0[mNCandidates] = trc0;
    mCandTr1[mNCandidates] = trc1;
    trc0Aux.glo2loc(mCrossings.xDCA[ic], mCrossings.yDCA[ic], xl,yl);
    if (!mCandTr0[mNCandidates].PropagateTo(xl,mBz)) {
      continue;
    }
    trc1Aux.glo2loc(mCrossings.xDCA[ic], mCrossings.yDCA[ic], xl,yl);
    if (!mCandTr1[mNCandidates].PropagateTo(xl,mBz)) {
      continue;
    }
    if (mUseAbsDCA ? processCandidateDCA(trc0Aux,trc1Aux) : processCandidateChi2(trc0Aux, trc1Aux)) {      
      mNCandidates++; // candidate validated
    }
  }
  return mNCandidates;
}

bool DCAFitter::processCandidateChi2(const TrackAuxPar& trc0Aux, const TrackAuxPar& trc1Aux)
{
  // find best chi2 (weighted DCA) of 2 tracks already propagated to their approximate vicinity
  mChi2[mNCandidates] = kDummyVal;
  AliExternalTrackParam& trc0 = mCandTr0[mNCandidates], &trc1 = mCandTr1[mNCandidates];
  TrackPointCovI& trcEI0 = mTrcEI0[mNCandidates], &trcEI1 =  mTrcEI1[mNCandidates];
  // get error matrices at initial point
  trcEI0.set(trc0);
  trcEI1.set(trc1);
  TrackCoefVtx& trCFVT0 = mTrCFVT0[mNCandidates], &trCFVT1 = mTrCFVT1[mNCandidates];  // get coefficients of PCA vs track points
  calcPCACoefs(trc0Aux, trcEI0, trc1Aux, trcEI1, trCFVT0, trCFVT1);
  TrackPoint tPnt0, tPnt1;                 // tracks local positions
  TrackPointDeriv2 tDer0, tDer1;           // their derivatives over track param X
  Derivatives deriv;                       // chi2 1st and 2nd derivatives over tracks params

  double chi2=0, chi2Prev = 0;
  int iter = 0;
  do {
    tPnt0.set(trc0); // tracks local positions  
    tPnt1.set(trc1); // at this iteration
    tDer0.set(trc0, mBz);
    tDer1.set(trc1, mBz);   // their derivatives over track param X    
    chi2Deriv(tPnt0, tDer0, trc0Aux, trcEI0, trCFVT0, tPnt1, tDer1, trc1Aux, trcEI1, trCFVT1, deriv);

    // do Newton-Rapson iteration with corrections = - dchi2/d{x0,x1} * [ d^2chi2/d{x0,x1}^2 ]^-1
    double detDer2 = deriv.dChidx0dx0*deriv.dChidx1dx1 - deriv.dChidx0dx1*deriv.dChidx0dx1;
    double detDer2I = 1./detDer2;
    double dX0 = -(deriv.dChidx0 * deriv.dChidx1dx1 - deriv.dChidx1 * deriv.dChidx0dx1)*detDer2I;
    double dX1 = -(deriv.dChidx1 * deriv.dChidx0dx0 - deriv.dChidx0 * deriv.dChidx0dx1)*detDer2I;
    trc0.PropagateParamOnlyTo(trc0.GetX() + dX0, mBz);
    trc1.PropagateParamOnlyTo(trc1.GetX() + dX1, mBz);
    
    double x,y,z;
    calcPCA(tPnt0,trCFVT0, tPnt1,trCFVT1, x,y,z);
    // if there are 2 XY crossings, fit with both of them used as a starting point may converge to
    // the same point. To stop asap this redundant step we abandon fit if it appears to be closer
    // to alternative crossing
    if (closerToAlternative(x,y)) {
      return false;
    } 
   
    chi2 = calcChi2(x,y,z, tPnt0, trc0Aux, trcEI0,  tPnt1, trc1Aux, trcEI1);

    if ( (TMath::Abs(dX0)<mMinParamChange && TMath::Abs(dX1)<mMinParamChange) ||
	 (iter && chi2/chi2Prev>mMinRelChi2Change) ) {
      break; // converged
    }
    chi2Prev = chi2;
  } while (++iter<mMaxIter);
  //
  if (chi2<mMaxChi2) {
    auto& pca = mPCA[mNCandidates];
    tPnt0.set(trc0); // tracks local positions  
    tPnt1.set(trc1); // 
    calcPCA(tPnt0, trCFVT0, tPnt1, trCFVT1, pca.x, pca.y, pca.z );
    mChi2[mNCandidates] = calcChi2(pca.x, pca.y, pca.z, tPnt0, trc0Aux, trcEI0,  tPnt1, trc1Aux, trcEI1);
    return true;
  }
  return false;
}

bool DCAFitter::processCandidateDCA(const TrackAuxPar& trc0Aux, const TrackAuxPar& trc1Aux)
{
  // find DCA of 2 tracks already propagated to its approximate vicinity W/O APPLYING ANY ERRORS
  // i.e. the absolute distance is minimized
  AliExternalTrackParam& trc0 = mCandTr0[mNCandidates], &trc1 = mCandTr1[mNCandidates];

  TrackPoint tPnt0, tPnt1;                 // tracks local positions
  TrackPointDeriv2 tDer0,tDer1;            // their derivatives over track param X
  Derivatives deriv;                       // chi2 1st and 2nd derivatives over tracks params

  double chi2=0, chi2Prev = 0;
  int iter = 0;
  do {
    tPnt0.set(trc0); // tracks local positions  
    tPnt1.set(trc1); // at this iteration
    tDer0.set(trc0, mBz);
    tDer1.set(trc1, mBz); // their derivatives over track param X
    
    DCADeriv(tPnt0, tDer0, trc0Aux, tPnt1, tDer1, trc1Aux, deriv);

    // do Newton-Rapson iteration with corrections = - dchi2/d{x0,x1} * [ d^2chi2/d{x0,x1}^2 ]^-1
    double detDer2 = deriv.dChidx0dx0*deriv.dChidx1dx1 - deriv.dChidx0dx1*deriv.dChidx0dx1;
    double detDer2I = 1./detDer2;    
    double dX0 = -(deriv.dChidx0 * deriv.dChidx1dx1 - deriv.dChidx1 * deriv.dChidx0dx1)*detDer2I;
    double dX1 = -(deriv.dChidx1 * deriv.dChidx0dx0 - deriv.dChidx0 * deriv.dChidx0dx1)*detDer2I;
    trc0.PropagateParamOnlyTo(trc0.GetX() + dX0, mBz);
    trc1.PropagateParamOnlyTo(trc1.GetX() + dX1, mBz);

    double x,y,z;
    tPnt0.set(trc0); // tracks local positions  
    tPnt1.set(trc1); // at this iteration
    calcPCA(tPnt0, trc0Aux, tPnt1, trc1Aux, x,y,z);
    // if there are 2 XY crossings, fit with both of them used as a starting point may converge to
    // the same point. To stop asap this redundant step we abandon fit if it appears to be closer
    // to alternative crossing
    if (closerToAlternative(x,y)) {
      return false;
    }

    chi2 = calcDCA(tPnt0, trc0Aux, tPnt1, trc1Aux);

    if ( (TMath::Abs(dX0)<mMinParamChange && TMath::Abs(dX1)<mMinParamChange) ||
	 (iter && chi2/chi2Prev>mMinRelChi2Change) ) {
      break; // converged
    }
    chi2Prev = chi2;
  } while (++iter<mMaxIter);
  //
  if (chi2<mMaxChi2) {
    auto& pca = mPCA[mNCandidates];
    tPnt0.set(trc0); // tracks local positions  
    tPnt1.set(trc1); //
    calcPCA(tPnt0, trc0Aux, tPnt1, trc1Aux, pca.x, pca.y, pca.z );
    mChi2[mNCandidates] = calcDCA(tPnt0, trc0Aux, tPnt1, trc1Aux);
    return true;
  }
  return false;
}


void DCAFitter::calcPCA(const AliExternalTrackParam& trc0, const TrackCoefVtx& trCFVT0,
			const AliExternalTrackParam& trc1, const TrackCoefVtx& trCFVT1,
			double &Vx, double &Vy, double &Vz) const
{
  // calculate the PCA (Vx,Vy,Vz) to 2 points in lab frame represented via local points coordinates
  Vx = trCFVT0.mXX*trc0.GetX() + trCFVT0.mXY*trc0.GetY() + trCFVT0.mXZ*trc0.GetX()
    +  trCFVT1.mXX*trc1.GetX() + trCFVT1.mXY*trc1.GetY() + trCFVT1.mXZ*trc1.GetX();
  Vy = trCFVT0.mYX*trc0.GetX() + trCFVT0.mYY*trc0.GetY() + trCFVT0.mYZ*trc0.GetX()
    +  trCFVT1.mYX*trc1.GetX() + trCFVT1.mYY*trc1.GetY() + trCFVT1.mYZ*trc1.GetX();
  Vz = trCFVT0.mZX*trc0.GetX() + trCFVT0.mZY*trc0.GetY() + trCFVT0.mZZ*trc0.GetX()
    +  trCFVT1.mZX*trc1.GetX() + trCFVT1.mZY*trc1.GetY() + trCFVT1.mZZ*trc1.GetX();
}

void DCAFitter::calcPCA(const TrackPoint& tPnt0, const TrackCoefVtx& trCFVT0,
			const TrackPoint& tPnt1, const TrackCoefVtx& trCFVT1,
			double &Vx, double &Vy, double &Vz) const
{
  // calculate the PCA (Vx,Vy,Vz) to 2 points in lab frame represented via local points coordinates
  Vx = trCFVT0.mXX*tPnt0.x + trCFVT0.mXY*tPnt0.y + trCFVT0.mXZ*tPnt0.z
    +  trCFVT1.mXX*tPnt1.x + trCFVT1.mXY*tPnt1.y + trCFVT1.mXZ*tPnt1.z;
  Vy = trCFVT0.mYX*tPnt0.x + trCFVT0.mYY*tPnt0.y + trCFVT0.mYZ*tPnt0.z
    +  trCFVT1.mYX*tPnt1.x + trCFVT1.mYY*tPnt1.y + trCFVT1.mYZ*tPnt1.z;
  Vz = trCFVT0.mZX*tPnt0.x + trCFVT0.mZY*tPnt0.y + trCFVT0.mZZ*tPnt0.z
    +  trCFVT1.mZX*tPnt1.x + trCFVT1.mZY*tPnt1.y + trCFVT1.mZZ*tPnt1.z;
}

void DCAFitter::calcPCA(const TrackPoint& tPnt0, const TrackAuxPar& trc0Aux, 
			const TrackPoint& tPnt1, const TrackAuxPar& trc1Aux, 
			double &Vx, double &Vy, double &Vz) const
{
  // calculate the PCA (Vx,Vy,Vz) to 2 points in lab frame represented via local points coordinates
  // w/o accounting for the errors of the points
  double xg,yg;
  trc0Aux.loc2glo(tPnt0.x, tPnt0.y, xg, yg);
  trc1Aux.loc2glo(tPnt1.x, tPnt1.y, Vx, Vy);
  Vx += xg;
  Vy += yg;
  Vz = 0.5*(tPnt0.z+tPnt1.z);
  Vx *= 0.5;
  Vy *= 0.5;
}

void DCAFitter::calcPCA(const AliExternalTrackParam& trc0, const TrackAuxPar& trc0Aux, 
			const AliExternalTrackParam& trc1, const TrackAuxPar& trc1Aux, 
			double &Vx, double &Vy, double &Vz) const
{
  // calculate the PCA (Vx,Vy,Vz) to 2 points in lab frame represented via local points coordinates
  // w/o accounting for the errors of the points
  double xg,yg;
  trc0Aux.loc2glo(trc0.GetX(), trc0.GetY(), xg, yg);
  trc1Aux.loc2glo(trc1.GetX(), trc1.GetY(), Vx, Vy);
  Vx += xg;
  Vy += yg;
  Vz = 0.5*(trc0.GetZ()+trc1.GetZ());
  Vx *= 0.5;
  Vy *= 0.5;
}

void DCAFitter::chi2Deriv(const TrackPoint& tPnt0, const TrackPointDeriv2& tDer0, const TrackAuxPar& trc0Aux, const TrackPointCovI& trcEI0, const TrackCoefVtx& trCFVT0,
			  const TrackPoint& tPnt1, const TrackPointDeriv2& tDer1, const TrackAuxPar& trc1Aux, const TrackPointCovI& trcEI1, const TrackCoefVtx& trCFVT1,
			  Derivatives& deriv) const
{
  // calculate 1st and 2nd derivatives of wighted DCA (chi2) over track parameters X
  TrackPoint vtx;
  calcPCA(tPnt0, trCFVT0, tPnt1, trCFVT1, vtx.x, vtx.y, vtx.z); // calculate PCA for current track-points positions
  // calculate residuals
  auto res0 = calcResid(tPnt0, trc0Aux, vtx);
  auto res1 = calcResid(tPnt1, trc1Aux, vtx);
  
  // res0.x,dy0,dz0 = x0 - Xl0, y0 - Yl0, ... , whith x0,y0,z0: track coords in its frame, and Xl0, Yl0, Zl0 : vertex rotated to same frame
  
  // aux params to minimize multiplications  
  double dx0s = res0.x*trcEI0.sxx;
  double dy0sz0t = res0.y*trcEI0.syy + res0.z*trcEI0.syz;
  double dy0tz0s = res0.y*trcEI0.syz + res0.z*trcEI0.szz;

  double dx1s = res1.x*trcEI1.sxx;
  double dy1sz1t = res1.y*trcEI1.syy + res1.z*trcEI1.syz;
  double dy1tz1s = res1.y*trcEI1.syz + res1.z*trcEI1.szz;

  //-------------------------------
  double xx0DtXYXZx0 = trCFVT0.mXX + trCFVT0.mXY*tDer0.dydx + trCFVT0.mXZ*tDer0.dzdx;
  double xx1DtXYXZx1 = trCFVT1.mXX + trCFVT1.mXY*tDer1.dydx + trCFVT1.mXZ*tDer1.dzdx;
  double yz0DtYYYZx0 = trCFVT0.mYX + trCFVT0.mYY*tDer0.dydx + trCFVT0.mYZ*tDer0.dzdx;
  double yz1DtYYYZx1 = trCFVT1.mYX + trCFVT1.mYY*tDer1.dydx + trCFVT1.mYZ*tDer1.dzdx;

  double DtXYXZx02 = trCFVT0.mXY*tDer0.d2ydx2 + trCFVT0.mXZ*tDer0.d2zdx2;
  double DtYYYZx02 = trCFVT0.mYY*tDer0.d2ydx2 + trCFVT0.mYZ*tDer0.d2zdx2;
  double DtXYXZx12 = trCFVT1.mXY*tDer1.d2ydx2 + trCFVT1.mXZ*tDer1.d2zdx2;
  double DtYYYZx12 = trCFVT1.mYY*tDer1.d2ydx2 + trCFVT1.mYZ*tDer1.d2zdx2;
    
  double FDdx0Dx0 = 1. - (trc0Aux.c*xx0DtXYXZx0 + trc0Aux.s*yz0DtYYYZx0);
  double FDdy0Dx0 = trc0Aux.s*xx0DtXYXZx0 - trc0Aux.c*yz0DtYYYZx0 + tDer0.dydx;
  double FDdz0Dx0 = -trCFVT0.mZX - trCFVT0.mZY*tDer0.dydx + tDer0.dzdx*(1. - trCFVT0.mZZ);
  
  double FDdx1Dx1 = 1. - (trc1Aux.c*xx1DtXYXZx1 + trc1Aux.s*yz1DtYYYZx1);
  double FDdy1Dx1 = trc1Aux.s*xx1DtXYXZx1 - trc1Aux.c*yz1DtYYYZx1 + tDer1.dydx;
  double FDdz1Dx1 = -trCFVT1.mZX - trCFVT1.mZY*tDer1.dydx + tDer1.dzdx*(1. - trCFVT1.mZZ);

  double FDdx0Dx1 = -(trc0Aux.c*xx1DtXYXZx1 + trc0Aux.s*yz1DtYYYZx1);
  double FDdy0Dx1 = trc0Aux.s*xx1DtXYXZx1 - trc0Aux.c*yz1DtYYYZx1;
  double FDdz0Dx1 = FDdz1Dx1 - tDer1.dzdx; // -trCFVT1.mZX - trCFVT1.mZY*tDer1.dydx - trCFVT1.mZZ*tDer1.dzdx;
    
  double FDdx1Dx0 = -(trc1Aux.c*xx0DtXYXZx0 + trc1Aux.s*yz0DtYYYZx0);
  double FDdy1Dx0 = trc1Aux.s*xx0DtXYXZx0 - trc1Aux.c*yz0DtYYYZx0;
  double FDdz1Dx0 = FDdz0Dx0 - tDer0.dzdx; //  -trCFVT0.mZX - trCFVT0.mZY*tDer0.dydx - trCFVT0.mZZ*tDer0.dzdx;

  double FDdx0Dx0x0 = -(trc0Aux.c*DtXYXZx02 + trc0Aux.s*DtYYYZx02);
  double FDdy0Dx0x0 = tDer0.d2ydx2 + trc0Aux.s*DtXYXZx02 - trc0Aux.c*DtYYYZx02;
  double FDdz0Dx0x0 = -trCFVT0.mZY*tDer0.d2ydx2 + tDer0.d2zdx2*(1. - trCFVT0.mZZ);

  double FDdx1Dx1x1 = -(trc1Aux.c*DtXYXZx12 + trc1Aux.s*DtYYYZx12);
  double FDdy1Dx1x1 = tDer1.d2ydx2 + trc1Aux.s*DtXYXZx12 - trc1Aux.c*DtYYYZx12;
  double FDdz1Dx1x1 = -trCFVT1.mZY*tDer1.d2ydx2 + tDer1.d2zdx2*(1. - trCFVT1.mZZ);

  double FDdx0Dx1x1 = -(trc0Aux.c*DtXYXZx12 + trc0Aux.s*DtYYYZx12);
  double FDdx1Dx0x0 = -(trc1Aux.c*DtXYXZx02 + trc1Aux.s*DtYYYZx02);

  double FDdy1Dx0x0 = trc1Aux.s*DtXYXZx02 - trc1Aux.c*DtYYYZx02;
  double FDdy0Dx1x1 = trc0Aux.s*DtXYXZx12 - trc0Aux.c*DtYYYZx12;

  double FDdz0Dx1x1 = FDdz1Dx1x1 - tDer1.d2zdx2; //  -(trCFVT1.mZY*tDer1.d2ydx2 + trCFVT1.mZZ*tDer1.d2zdx2);
  double FDdz1Dx0x0 = FDdz0Dx0x0 - tDer0.d2zdx2; //  -(trCFVT0.mZY*tDer0.d2ydx2 + trCFVT0.mZZ*tDer0.d2zdx2);
  
  double FD00YYYZ = FDdy0Dx0*trcEI0.syy + FDdz0Dx0*trcEI0.syz;
  double FD11YYYZ = FDdy1Dx1*trcEI1.syy + FDdz1Dx1*trcEI1.syz;
  double FD10YYYZ = FDdy1Dx0*trcEI1.syy + FDdz1Dx0*trcEI1.syz;
  double FD01YYYZ = FDdy0Dx1*trcEI0.syy + FDdz0Dx1*trcEI0.syz;
  double FD00YZZZ = FDdy0Dx0*trcEI0.syz + FDdz0Dx0*trcEI0.szz;
  double FD11YZZZ = FDdy1Dx1*trcEI1.syz + FDdz1Dx1*trcEI1.szz;

  double FD10YZZZ = FDdy1Dx0*trcEI1.syz + FDdz1Dx0*trcEI1.szz;
  double FD01YZZZ = FDdy0Dx1*trcEI0.syz + FDdz0Dx1*trcEI0.szz;
  
  // 1st derivatives over track params x
  deriv.dChidx0 = dx0s*FDdx0Dx0 + dx1s*FDdx1Dx0 + dy0sz0t*FDdy0Dx0 + dy1sz1t*FDdy1Dx0 + dy0tz0s*FDdz0Dx0 + dy1tz1s*FDdz1Dx0;
  deriv.dChidx1 = dx1s*FDdx1Dx1 + dx0s*FDdx0Dx1 + dy1sz1t*FDdy1Dx1 + dy0sz0t*FDdy0Dx1 + dy1tz1s*FDdz1Dx1 + dy0tz0s*FDdz0Dx1;
  
  
  // 2nd derivative over track params x
  deriv.dChidx0dx0 = dx0s*FDdx0Dx0x0 + dx1s*FDdx1Dx0x0 + FDdy0Dx0x0*dy0sz0t + FDdy1Dx0x0*dy1sz1t + FDdz0Dx0x0*dy0tz0s + FDdz1Dx0x0*dy1tz1s +
    FDdx0Dx0*FDdx0Dx0*trcEI0.sxx + FDdx1Dx0*FDdx1Dx0*trcEI1.sxx +
    FDdy0Dx0*FD00YYYZ + FDdy1Dx0*FD10YYYZ + FDdz0Dx0*FD00YZZZ + FDdz1Dx0*FD10YZZZ;
  
  deriv.dChidx1dx1 = dx1s*FDdx1Dx1x1 + dx0s*FDdx0Dx1x1 + FDdy1Dx1x1*dy1sz1t + FDdy0Dx1x1*dy0sz0t + FDdz1Dx1x1*dy1tz1s + FDdz0Dx1x1*dy0tz0s +
    FDdx1Dx1*FDdx1Dx1*trcEI1.sxx + FDdx0Dx1*FDdx0Dx1*trcEI0.sxx +
    FDdy1Dx1*FD11YYYZ + FDdy0Dx1*FD01YYYZ + FDdz1Dx1*FD11YZZZ + FDdz0Dx1*FD01YZZZ;

  // N.B.: cross-derivatice
  //  FDdx0Dx0x1 -> 0, FDdy0Dx0x1 -> 0, FDdz0Dx0x1 -> 0, FDdx1Dx1x0 -> 0;
  //  FDdy1Dx1x0 -> 0, FDdz1Dx1x0 -> 0, FDdx0Dx1x0 -> 0, FDdy0Dx1x0 -> 0;
  //  FDdz0Dx1x0 -> 0, FDdx1Dx0x1 -> 0, FDdy1Dx0x1 -> 0, FDdz1Dx0x1 -> 0;
  
  deriv.dChidx0dx1 =
    //    this part is = 0 due to the N.B.
    //    dx0s*FDdx0Dx0x1 + dx1s*FDdx1Dx0x1 + dy0s*FDdy0Dx0x1 + dy1s*FDdy1Dx0x1 + dz0s*FDdz0Dx0x1 + dz1s*FDdz1Dx0x1 + 
    //    dy0t*FDdz0Dx0x1 + dy1t*FDdz1Dx0x1 +dz0t*FDdy0Dx0x1 + dz1t*FDdy1Dx0x1 +     
    FDdx0Dx0*FDdx0Dx1*trcEI0.sxx +  FDdx1Dx0*FDdx1Dx1*trcEI1.sxx +
    FDdy0Dx0*(FDdy0Dx1*trcEI0.syy + FDdz0Dx1*trcEI0.syz) +
    FDdy1Dx0*(FDdy1Dx1*trcEI1.syy + FDdz1Dx1*trcEI1.syz) +    
    FDdz0Dx0*(FDdy0Dx1*trcEI0.syz + FDdz0Dx1*trcEI0.szz) +
    FDdz1Dx0*(FDdy1Dx1*trcEI1.syz + FDdz1Dx1*trcEI1.szz);
  //
}

void DCAFitter::DCADeriv(const TrackPoint& tPnt0, const TrackPointDeriv2& tDer0, const TrackAuxPar& trc0Aux, 
			 const TrackPoint& tPnt1, const TrackPointDeriv2& tDer1, const TrackAuxPar& trc1Aux, 
			 Derivatives& deriv) const
{
  // DCA derivative calculation with Chi2 defined as the absolule distance (no errors applied)
  double cosDA = trc0Aux.c*trc1Aux.c + trc0Aux.s*trc1Aux.s; // cos(A0-A1)
  double sinDA = trc0Aux.s*trc1Aux.c - trc0Aux.c*trc1Aux.s; // sin(A0-A1)  
  double dx01 = tPnt0.x - (tPnt1.x*cosDA + tPnt1.y*sinDA);
  double dy01 = tPnt0.y - (tPnt1.y*cosDA - tPnt1.x*sinDA);
  double dx10 = tPnt1.x - (tPnt0.x*cosDA - tPnt0.y*sinDA);
  double dy10 = tPnt1.y - (tPnt0.y*cosDA + tPnt0.x*sinDA);
  double dz = tPnt0.z - tPnt1.z;

  // 1st derivatives over track params x
  deriv.dChidx0 = 0.5*(tDer0.dydx*dy01 + dx01 + dz*tDer0.dzdx);
  deriv.dChidx1 = 0.5*(tDer1.dydx*dy10 + dx10 - dz*tDer1.dzdx);

  // 2nd derivative over track params x
  deriv.dChidx0dx0 = 0.5*(1. + tDer0.dydx*tDer0.dydx + tDer0.dzdx*tDer0.dzdx + dz*tDer0.d2zdx2 + tDer0.d2ydx2*dy01);
  deriv.dChidx1dx1 = 0.5*(1. + tDer1.dydx*tDer1.dydx + tDer1.dzdx*tDer1.dzdx - dz*tDer1.d2zdx2 + tDer1.d2ydx2*dy10);
  deriv.dChidx0dx1 = 0.5*(sinDA*(tDer0.dydx - tDer1.dydx) - cosDA*(1. + tDer0.dydx*tDer1.dydx) - tDer0.dzdx*tDer1.dzdx);
}


void DCAFitter::calcPCACoefs(const TrackAuxPar& trc0Aux, const TrackPointCovI& trcEI0,
			     const TrackAuxPar& trc1Aux, const TrackPointCovI& trcEI1,
			     TrackCoefVtx& trCFVT0, TrackCoefVtx& trCFVT1) const
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
  double axx = trc0Aux.cc*trcEI0.sxx + trc1Aux.cc*trcEI1.sxx + trc0Aux.ss*trcEI0.syy + trc1Aux.ss*trcEI1.syy;
  double axy = trc0Aux.cs*(trcEI0.sxx - trcEI0.syy) + trc1Aux.cs*(trcEI1.sxx - trcEI1.syy);
  double axz = -(trc0Aux.s*trcEI0.syz + trc1Aux.s*trcEI1.syz);
  double ayy = trc0Aux.ss*trcEI0.sxx + trc1Aux.ss*trcEI1.sxx + trc0Aux.cc*trcEI0.syy + trc1Aux.cc*trcEI1.syy; // = (trcEI0.sxx + trcEI1.sxx + trcEI0.syy + trcEI1.syy) - axx
  double ayz = trc0Aux.c*trcEI0.syz + trc1Aux.c*trcEI1.syz;
  double azz = trcEI0.szz + trcEI1.szz;
  //
  // define some aux variables
  double axxyy = axx*ayy, axxzz = axx*azz, axxyz = axx*ayz, 
    axyxy = axy*axy, axyxz = axy*axz, axyyz = axy*ayz, axyzz = axy*azz, 
    axzxz = axz*axz, axzyy = axz*ayy, axzyz = axz*ayz,
    ayzyz = ayz*ayz,  ayyzz = ayy*azz;
  double dAxxyyAxyxy = axxyy - axyxy, dAxyyzAxzyy = axyyz - axzyy, dAxyxzAxxyz = axyxz - axxyz;
  double dAxzyzAxyzz = axzyz - axyzz, dAyyzzAyzyz = ayyzz - ayzyz, dAxxzzAxzxz = axxzz - axzxz;
  double det = -dAxyyzAxzyy*axz + dAxyxzAxxyz*ayz + dAxxyyAxyxy*azz, detI = 1./det;
  
  double dfxPCS0 = dAyyzzAyzyz*trc0Aux.c + dAxzyzAxyzz*trc0Aux.s, dfxQCS0 = dAxzyzAxyzz*trc0Aux.c - dAyyzzAyzyz*trc0Aux.s;
  double dfyPCS0 = dAxzyzAxyzz*trc0Aux.c + dAxxzzAxzxz*trc0Aux.s, dfyQCS0 = dAxxzzAxzxz*trc0Aux.c - dAxzyzAxyzz*trc0Aux.s;
  double dfzPCS0 = dAxyyzAxzyy*trc0Aux.c + dAxyxzAxxyz*trc0Aux.s, dfzQCS0 = dAxyxzAxxyz*trc0Aux.c - dAxyyzAxzyy*trc0Aux.s;

  double dfxPCS1 = dAyyzzAyzyz*trc1Aux.c + dAxzyzAxyzz*trc1Aux.s, dfxQCS1 = dAxzyzAxyzz*trc1Aux.c - dAyyzzAyzyz*trc1Aux.s;
  double dfyPCS1 = dAxzyzAxyzz*trc1Aux.c + dAxxzzAxzxz*trc1Aux.s, dfyQCS1 = dAxxzzAxzxz*trc1Aux.c - dAxzyzAxyzz*trc1Aux.s;
  double dfzPCS1 = dAxyyzAxzyy*trc1Aux.c + dAxyxzAxxyz*trc1Aux.s, dfzQCS1 = dAxyxzAxxyz*trc1Aux.c - dAxyyzAxzyy*trc1Aux.s;
  //
  trCFVT0.mXX = detI*(dfxPCS0*trcEI0.sxx);
  trCFVT0.mXY = detI*(dfxQCS0*trcEI0.syy + dAxyyzAxzyy*trcEI0.syz);
  trCFVT0.mXZ = detI*(dfxQCS0*trcEI0.syz + dAxyyzAxzyy*trcEI0.szz);
  
  trCFVT0.mYX = detI*(dfyPCS0*trcEI0.sxx);
  trCFVT0.mYY = detI*(dfyQCS0*trcEI0.syy + dAxyxzAxxyz*trcEI0.syz);
  trCFVT0.mYZ = detI*(dfyQCS0*trcEI0.syz + dAxyxzAxxyz*trcEI0.szz);
  
  trCFVT0.mZX = detI*(dfzPCS0*trcEI0.sxx);
  trCFVT0.mZY = detI*(dfzQCS0*trcEI0.syy + dAxxyyAxyxy*trcEI0.syz);
  trCFVT0.mZZ = detI*(dfzQCS0*trcEI0.syz + dAxxyyAxyxy*trcEI0.szz);
  
  trCFVT1.mXX = detI*(dfxPCS1*trcEI1.sxx);
  trCFVT1.mXY = detI*(dfxQCS1*trcEI1.syy + dAxyyzAxzyy*trcEI1.syz);
  trCFVT1.mXZ = detI*(dfxQCS1*trcEI1.syz + dAxyyzAxzyy*trcEI1.szz);
  
  trCFVT1.mYX = detI*(dfyPCS1*trcEI1.sxx);
  trCFVT1.mYY = detI*(dfyQCS1*trcEI1.syy + dAxyxzAxxyz*trcEI1.syz);
  trCFVT1.mYZ = detI*(dfyQCS1*trcEI1.syz + dAxyxzAxxyz*trcEI1.szz);
  
  trCFVT1.mZX = detI*(dfzPCS1*trcEI1.sxx);
  trCFVT1.mZY = detI*(dfzQCS1*trcEI1.syy + dAxxyyAxyxy*trcEI1.syz);
  trCFVT1.mZZ = detI*(dfzQCS1*trcEI1.syz + dAxxyyAxyxy*trcEI1.szz);
  
}

double DCAFitter::calcChi2(double Vx, double Vy, double Vz,
			   const TrackPoint& tPnt0, const TrackAuxPar& trc0Aux, const TrackPointCovI& trcEI0,
			   const TrackPoint& tPnt1, const TrackAuxPar& trc1Aux, const TrackPointCovI& trcEI1) const
{
  double chi2 = 0;
  double xl,yl,dy,dz;
  trc0Aux.glo2loc(Vx,Vy,xl,yl);
  dy = tPnt0.y - yl;
  dz = tPnt0.z - Vz;
  chi2 += dy*dy*trcEI0.syy + dz*dz*trcEI0.szz + 2.*dy*dz*trcEI0.syz;
  trc1Aux.glo2loc(Vx,Vy,xl,yl);
  dy = tPnt1.y - yl;
  dz = tPnt1.z - Vz;
  chi2 += dy*dy*trcEI1.syy + dz*dz*trcEI1.szz + 2.*dy*dz*trcEI1.syz;
  return chi2;
}


double DCAFitter::calcDCA(double Vx, double Vy, double Vz,
			  const TrackPoint& tPnt0, const TrackAuxPar& trc0Aux, const TrackPoint& tPnt1, const TrackAuxPar& trc1Aux) const
{
  // calculate distance (non-weighted) of closest approach of 2 points in their local frame
  // (long way, see alternative getDCA w/o explicit vertex calculation)
  double chi2 = 0;
  double xl,yl,dy,dz;
  trc0Aux.glo2loc(Vx,Vy,xl,yl);
  dy = tPnt0.y - yl;
  dz = tPnt0.z - Vz;
  chi2 += dy*dy + dz*dz;
  trc1Aux.glo2loc(Vx,Vy,xl,yl);
  dy = tPnt1.y - yl;
  dz = tPnt1.z - Vz;
  chi2 += dy*dy + dz*dz;
  return chi2;
}


double DCAFitter::calcDCA( const TrackPoint& tPnt0, const TrackAuxPar& trc0Aux, const TrackPoint& tPnt1, const TrackAuxPar& trc1Aux) const
{
  // calculate distance (non-weighted) of closest approach of 2 points in their local frame
  double chi2 = 0;
  double cosDA = trc0Aux.c*trc1Aux.c + trc0Aux.s*trc1Aux.s; // cos(A0-A1)
  double sinDA = trc0Aux.s*trc1Aux.c - trc0Aux.c*trc1Aux.s; // sin(A0-A1)
  double dx = tPnt0.x - tPnt1.x, dy = tPnt0.y - tPnt1.y, dz = tPnt0.z - tPnt1.z; 
  chi2 = 0.5*(dx*dx + dy*dy + dz*dz) + (1.-cosDA)*(tPnt0.x*tPnt1.x + tPnt0.y*tPnt1.y) + sinDA*(tPnt0.y*tPnt1.x-tPnt1.y*tPnt0.x);
  
  return chi2;
}

bool DCAFitter::closerToAlternative(float x, float y) const
{
  // check if the point x,y is closer to the seeding XY point being tested or to alternative see (if any)
  if (mCrossIDAlt<0) {
    return false;
  }
  float dxCur = x - mCrossings.xDCA[mCrossITCur], dyCur = y - mCrossings.yDCA[mCrossITCur];
  float dxAlt = x - mCrossings.xDCA[mCrossIDAlt], dyAlt = y - mCrossings.yDCA[mCrossIDAlt];
  return dxCur*dxCur + dyCur*dyCur > dxAlt*dxAlt + dyAlt*dyAlt;
}
