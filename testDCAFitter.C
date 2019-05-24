
void testDCAFitter()
{
  double bz = 5.0;

  // create V0 
  AliExternalTrackParam t0;
  t0.GetParameter()[4] = 1.;
  t0.GetParameter()[3] = 0.6;
  t0.GetParameter()[2] = -0.2;
  double c[15]={1e-6,0,1e-6,0,0,1e-6,0,0,0,1e-6,0,0,0,0,1e-5};
  memcpy(t0.GetCovariance(),c,15*sizeof(double));
  t0.PropagateTo(10,bz);

  double xyz[3];
  t0.GetXYZ(xyz);
  printf("true vertex : %+e %+e %+e\n",xyz[0],xyz[1],xyz[2]);

  AliExternalTrackParam tA(t0), tB(t0);
  tB.GetParameter()[4] = -tA.GetParameter()[4];
  tB.GetParameter()[3] = -0.5*tA.GetParameter()[3];
  tB.GetParameter()[2] = -0.2*tA.GetParameter()[2];

  // randomize the tracks positions
  tA.GetParameter()[0] += gRandom->Gaus()*TMath::Sqrt(tA.GetSigmaY2()); 
  tA.GetParameter()[1] += gRandom->Gaus()*TMath::Sqrt(tA.GetSigmaZ2()); 
  tB.GetParameter()[0] += gRandom->Gaus()*TMath::Sqrt(tB.GetSigmaY2()); 
  tB.GetParameter()[1] += gRandom->Gaus()*TMath::Sqrt(tB.GetSigmaZ2()); 

  tA.Rotate(-0.3);
  tB.Rotate(0.3);

  printf("True track params at PCA:\n");
  tA.Print();
  tB.Print();
  
  tA.PropagateTo(tA.GetX()+5.,bz);
  tB.PropagateTo(tB.GetX()+8.,bz);

  DCAFitter df(bz,10.);

  printf("\n\nTesting with abs DCA minimization\n");
  df.setUseAbsDCA(true);

  // we may supply track directly to the fitter (they are not modified)
  int nCand = df.process(tA,tB);
  // we can have up to 2 candidates
  for (int ic=0;ic<nCand;ic++) {
    const DCAFitter::TrackPoint& vtx = df.getPCACandidate(ic);
    const AliExternalTrackParam& trc0 = df.getTrack0(ic), &trc1 = df.getTrack1(ic); // track parameters at V0
    printf("Candidate %d: DCA:%+e Vtx: %+e %+e %+e [Diff to true: %+e %+e %+e]\n", ic, df.getChi2AtPCACandidate(ic),
	   vtx.x,vtx.y,vtx.z, vtx.x-xyz[0], vtx.y-xyz[1], vtx.z-xyz[2] );
    printf("Track X-parameters at PCA: %+e %+e\n",trc0.GetX(), trc1.GetX());
    trc0.Print();
    trc1.Print();
  } 
  
  printf("\n\nTesting with weighted DCA minimization\n");
  df.setUseAbsDCA(false);
  
  // we can also precalculate some track parameters: useful in a tight loops with the same track tested in many V0s
  DCAFitter::TrackAuxPar trc0Aux(tA,bz), trc1Aux(tB,bz);
  nCand = df.process(tA,tB);
  for (int ic=0;ic<nCand;ic++) {
    const DCAFitter::TrackPoint& vtx = df.getPCACandidate(ic);
    const AliExternalTrackParam& trc0 = df.getTrack0(ic), &trc1 = df.getTrack1(ic); // track parameters at V0
    printf("Candidate %d: DCA:%+e Vtx: %+e %+e %+e [Diff to true: %+e %+e %+e]\n", ic, df.getChi2AtPCACandidate(ic),
	   vtx.x,vtx.y,vtx.z, vtx.x-xyz[0], vtx.y-xyz[1], vtx.z-xyz[2] );
    printf("Track X-parameters at PCA: %+e %+e\n",trc0.GetX(), trc1.GetX());
    trc0.Print();
    trc1.Print();
  }

}
