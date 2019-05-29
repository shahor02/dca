#include "DCAFitter.h"
#include "TRandom.h"

void testO2_dcaFitter()
{
  double bz = 5.0;

  // create V0 
  std::array<float,15> cv={1e-6,0,1e-6,0,0,1e-6,0,0,0,1e-6,0,0,0,0,1e-5};
  std::array<float,5> pr={0.,0.,-0.2,0.6,1.};  
  Track t0(0.,0.,pr,cv);
  t0.propagateTo(10,bz);

  std::array<float,3> xyz;
  t0.getXYZGlo(xyz);
  printf("true vertex : %+e %+e %+e\n",xyz[0],xyz[1],xyz[2]);

  Track tA(t0), tB(t0);
  tB.setParam(-tA.getParam(4), 4);
  tB.setParam(-0.5*tA.getParam(3), 3);
  tB.setParam(-0.2*tA.getParam(2), 2);

  // randomize the tracks positions
  tA.setParam(tA.getParam(0) + gRandom->Gaus()*TMath::Sqrt(tA.getSigmaY2()), 0); 
  tA.setParam(tA.getParam(1) + gRandom->Gaus()*TMath::Sqrt(tA.getSigmaZ2()), 1); 
  tB.setParam(tB.getParam(0) + gRandom->Gaus()*TMath::Sqrt(tB.getSigmaY2()), 0); 
  tB.setParam(tB.getParam(1) + gRandom->Gaus()*TMath::Sqrt(tB.getSigmaZ2()), 1);
  
  tA.rotate(-0.3);
  tB.rotate(0.3);

  printf("True track params at PCA:\n");
  tA.print();
  tB.print();
  
  tA.propagateTo(tA.getX()+5.,bz);
  tB.propagateTo(tB.getX()+8.,bz);

  DCAFitter df(bz,10.);

  
  df.setUseAbsDCA(true);
  // we may supply track directly to the fitter (they are not modified)
  int nCand = df.process(tA,tB);
  printf("\n\nTesting with abs DCA minimization: %d candidates found\n",nCand);
  // we can have up to 2 candidates
  for (int ic=0;ic<nCand;ic++) {
    const DCAFitter::Triplet& vtx = df.getPCACandidate(ic);
    const Track& trc0 = df.getTrack0(ic), &trc1 = df.getTrack1(ic); // track parameters at V0
    printf("Candidate %d: DCA:%+e Vtx: %+e %+e %+e [Diff to true: %+e %+e %+e]\n", ic, df.getChi2AtPCACandidate(ic),
	   vtx.x,vtx.y,vtx.z, vtx.x-xyz[0], vtx.y-xyz[1], vtx.z-xyz[2] );
    printf("Track X-parameters at PCA: %+e %+e\n",trc0.getX(), trc1.getX());
    trc0.print();
    trc1.print();
  } 

  
  df.setUseAbsDCA(false);  
  // we can also precalculate some track parameters: useful in a tight loops with the same track tested in many V0s
  DCAFitter::TrcAuxPar trc0Aux(tA,bz), trc1Aux(tB,bz);
  nCand = df.process(tA,tB);
  printf("\n\nTesting with weighted DCA minimization: %d candidates found\n",nCand);
  for (int ic=0;ic<nCand;ic++) {
    const DCAFitter::Triplet& vtx = df.getPCACandidate(ic);
    const Track& trc0 = df.getTrack0(ic), &trc1 = df.getTrack1(ic); // track parameters at V0
    printf("Candidate %d: DCA:%+e Vtx: %+e %+e %+e [Diff to true: %+e %+e %+e]\n", ic, df.getChi2AtPCACandidate(ic),
	   vtx.x,vtx.y,vtx.z, vtx.x-xyz[0], vtx.y-xyz[1], vtx.z-xyz[2] );
    printf("Track X-parameters at PCA: %+e %+e\n",trc0.getX(), trc1.getX());
    trc0.print();
    trc1.print();
  }

}
