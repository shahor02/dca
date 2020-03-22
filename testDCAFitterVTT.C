#if !(defined(__CLING__)  || defined(__CINT__)) || defined(__ROOTCLING__) || defined(__ROOTCINT__)
#include "DCAFitter.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"
#include "TObjArray.h"
#endif


void testDCAFitterVTT()
{
  double bz = 5.0;

  // create V0
  double alp0 = 0.5;
  double x0 = 10.;
  double p[5] = {0,0,-0.2,0.6,1.0};
  double c[15]={1e-6,0,1e-6,0,0,1e-6,0,0,0,1e-6,0,0,0,0,1e-5};
  AliExternalTrackParam tA(x0,alp0, p, c);
  p[4] = -p[4];
  p[3] = -0.5*p[3];
  p[2] = -0.2*p[2];
  AliExternalTrackParam tB(x0,alp0, p, c);  
  p[4] = -p[4]/2;
  p[3] = 0.5*p[3];
  p[2] = 0.2*p[2];
  AliExternalTrackParam tC(x0,alp0, p, c);  

  double xyz[3];
  tA.GetXYZ(xyz);
  printf("true vertex : %+e %+e %+e\n",xyz[0],xyz[1],xyz[2]);

  tA.Rotate(alp0 - 0.3);
  tB.Rotate(alp0 + 0.3);
  tB.Rotate(alp0 - 0.1);

  printf("True track params at PCA:\n");
  tA.Print();
  tB.Print();
  tC.Print();
  
  tA.PropagateTo(tA.GetX()+5.,bz);
  tB.PropagateTo(tB.GetX()+8.,bz);
  tC.PropagateTo(tC.GetX()-2.,bz);

  DCAFitter df(bz,10.);

  AliVertexerTracks VTT(bz);
  AliESDVertex esv0;
  esv0.SetXv(xyz[0]);
  esv0.SetYv(xyz[1]);
  esv0.SetZv(xyz[2]);
  TObjArray trkArray;
  trkArray.Add(&tA);
  trkArray.Add(&tB);
  UShort_t idDumm[4]={0,1,2,3};
  
  df.setUseAbsDCA(true);
  // we may supply track directly to the fitter (they are not modified)
  double sm = 0;
  TStopwatch sw0;
  for (int i=0;i<1000000;i++) {
    sm += df.process(tA,tB);
  }
  sw0.Stop();
  printf("time for DCA2 wgh %e\n",sm);
  sw0.Print();

  sm = 0;
  sw0.Start();
  for (int i=0;i<1000000;i++) {
    VTT.SetVtxStart(&esv0);
    AliESDVertex *vertexESD = VTT.VertexForSelectedTracks(&trkArray,&idDumm[0]);
    sm += vertexESD->GetZ();
  }
  sw0.Stop();
  printf("time for VertTrack abs %e\n",sm);
  sw0.Print();

  
  trkArray.Add(&tC);
  sm = 0;
  sw0.Start();
  for (int i=0;i<1000000;i++) {
    VTT.SetVtxStart(&esv0);
    AliESDVertex *vertexESD = VTT.VertexForSelectedTracks(&trkArray,&idDumm[0]);
    sm += vertexESD->GetZ();
  }
  sw0.Stop();
  printf("time for 3Prong VertTrack abs %e\n",sm);
  sw0.Print();

  
  
  int nCand = df.process(tA,tB);
  printf("\n\nTesting with abs DCA minimization: %d candidates found\n",nCand);
  // we can have up to 2 candidates
  for (int ic=0;ic<nCand;ic++) {
    const DCAFitter::Triplet& vtx = df.getPCACandidate(ic);
    float dx = vtx.x - xyz[0], dy = vtx.y - xyz[1], dz = vtx.z - xyz[2];
    float dst = TMath::Sqrt(dx * dx + dy * dy + dz * dz);

    const AliExternalTrackParam& trc0 = df.getTrack0(ic), &trc1 = df.getTrack1(ic); // track parameters at V0
    printf("Candidate %d: DCA:%+e Vtx: %+e %+e %+e [Diff to true: %+e %+e %+e -> %+e]\n",
           ic, df.getChi2AtPCACandidate(ic), vtx.x, vtx.y, vtx.z, dx, dy, dz, dst);

    float dist2 = DCAFitter::getDistance2(trc0,trc1);
    float dist2v = DCAFitter::getDistance2(vtx.x, vtx.y, vtx.z, trc0,trc1);
    printf("Dist^2: track-track: %+.3e | track-v0-track: %+.3e\n",dist2, dist2v);
    
    printf("Track X-parameters at PCA: %+e %+e\n",trc0.GetX(), trc1.GetX());
    trc0.Print();
    trc1.Print();
  } 

  
  df.setUseAbsDCA(false);  
  // we can also precalculate some track parameters: useful in a tight loops with the same track tested in many V0s
  DCAFitter::TrcAuxPar trc0Aux(tA,bz), trc1Aux(tB,bz);

  sm = 0;
  sw0.Start();
  for (int i=0;i<1000000;i++) {
    sm += df.process(tA,tB);
  }
  sw0.Stop();
  printf("time for DCA2 wgh %e\n",sm);
  sw0.Print();


  nCand = df.process(tA,tB);
  printf("\n\nTesting with weighted DCA minimization: %d candidates found\n",nCand);
  for (int ic=0;ic<nCand;ic++) {
    const DCAFitter::Triplet& vtx = df.getPCACandidate(ic);
    float dx = vtx.x - xyz[0], dy = vtx.y - xyz[1], dz = vtx.z - xyz[2];
    float dst = TMath::Sqrt(dx * dx + dy * dy + dz * dz);
    const AliExternalTrackParam& trc0 = df.getTrack0(ic), &trc1 = df.getTrack1(ic); // track parameters at V0
    printf("Candidate %d: DCA:%+e Vtx: %+e %+e %+e [Diff to true: %+e %+e %+e -> %e]\n",
           ic, df.getChi2AtPCACandidate(ic), vtx.x, vtx.y, vtx.z, dx, dy, dz, dst);
    
    float dist2 = DCAFitter::getDistance2(trc0,trc1);
    float dist2v = DCAFitter::getDistance2(vtx.x, vtx.y, vtx.z, trc0,trc1);
    printf("Dist^2: track-track: %+.3e | track-v0-track: %+.3e\n",dist2, dist2v);

    trc0.Print();
    trc1.Print();
  }

}
