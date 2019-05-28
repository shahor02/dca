DCA of 2 AliExternalTrackParam tracks, some documentation files on
https://alice.its.cern.ch/jira/browse/O2-712

Supports both DCA (unweighted) and Chi2 (weighted DCA) minimization.
See testDCAFitter for example of usage:
aliroot
root [0] .L DCAFitter.cxx+
root [1] .x testDCAFitter.C

Starts by finding analytical XY crossings (at most 2) and using them as 
seeds for non-linear minimization (the track errors used for Chi2 are 
evaluated only once at these seeding points).
In case the PCA being fitted with one seed appears to be closer to 
alternative seed, the minimizaiton with this seed is abandoned.
In some cases 2 seeds may produce 2 distinct PCA candidates, therefore one should
check the DCAFitter::getNCandidates() result.

Files dcaMathFin.nb (dcaMathFinNoErr.nb) are Mathematica derivations
of the chi2 (or DCA) 1st and 2nd derivatives used in the code for
Newton-Rapson iterations.

----------------------------

To do: 
1. Port to o2::track model
2. At the moment the propogation is done the simple PropagateTo method using constant bz field.
Need to inteface the fitter to o2::base::Propagator and consider material corrections.
3. More elaborate analysis of seeding XY points is needed, in most of cases with 2 crossings
one of them can be rejected even w/o testing the DCA.
