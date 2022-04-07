\ Model tree-optimizer
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0.1666666666666667 Z[KSCYng] + 0.1666666666666667 Z[SNVAng]
   + 0.1666666666666667 Z[STTLng] + 0.1666666666666667 Z[ATLAng]
   + 0.1666666666666667 Z[CHINng] + 0.1666666666666667 Z[IPLSng]
   + 0.1666666666666667 Z[NYCMng] + 0.1666666666666667 Z[ATLAM5]
   + 0.1666666666666667 Z[LOSAng] + 0.1666666666666667 Z[DNVRng]
   + 0.1666666666666667 Z[WASHng] + 0.1666666666666667 Z[HSTNng]
   + 0.1666666666666667 Z[KSCYng] + 0.1666666666666667 Z[SNVAng]
   + 0.1666666666666667 Z[STTLng] + 0.1666666666666667 Z[ATLAng]
   + 0.1666666666666667 Z[CHINng] + 0.1666666666666667 Z[IPLSng]
   + 0.1666666666666667 Z[NYCMng] + 0.1666666666666667 Z[ATLAM5]
   + 0.1666666666666667 Z[LOSAng] + 0.1666666666666667 Z[DNVRng]
   + 0.1666666666666667 Z[WASHng] + 0.1666666666666667 Z[HSTNng]
   + 0.1666666666666667 Z[KSCYng] + 0.1666666666666667 Z[SNVAng]
   + 0.1666666666666667 Z[STTLng] + 0.1666666666666667 Z[ATLAng]
   + 0.1666666666666667 Z[CHINng] + 0.1666666666666667 Z[IPLSng]
   + 0.1666666666666667 Z[NYCMng] + 0.1666666666666667 Z[ATLAM5]
   + 0.1666666666666667 Z[LOSAng] + 0.1666666666666667 Z[DNVRng]
   + 0.1666666666666667 Z[WASHng] + 0.1666666666666667 Z[HSTNng]
   + 0.1666666666666667 Z[KSCYng] + 0.1666666666666667 Z[SNVAng]
   + 0.1666666666666667 Z[STTLng] + 0.1666666666666667 Z[ATLAng]
   + 0.1666666666666667 Z[CHINng] + 0.1666666666666667 Z[IPLSng]
   + 0.1666666666666667 Z[NYCMng] + 0.1666666666666667 Z[ATLAM5]
   + 0.1666666666666667 Z[LOSAng] + 0.1666666666666667 Z[DNVRng]
   + 0.1666666666666667 Z[WASHng] + 0.1666666666666667 Z[HSTNng]
   + 0.1666666666666667 Z[KSCYng] + 0.1666666666666667 Z[SNVAng]
   + 0.1666666666666667 Z[STTLng] + 0.1666666666666667 Z[ATLAng]
   + 0.1666666666666667 Z[CHINng] + 0.1666666666666667 Z[IPLSng]
   + 0.1666666666666667 Z[NYCMng] + 0.1666666666666667 Z[ATLAM5]
   + 0.1666666666666667 Z[LOSAng] + 0.1666666666666667 Z[DNVRng]
   + 0.1666666666666667 Z[WASHng] + 0.1666666666666667 Z[HSTNng]
   + 0.1666666666666667 Z[KSCYng] + 0.1666666666666667 Z[SNVAng]
   + 0.1666666666666667 Z[STTLng] + 0.1666666666666667 Z[ATLAng]
   + 0.1666666666666667 Z[CHINng] + 0.1666666666666667 Z[IPLSng]
   + 0.1666666666666667 Z[NYCMng] + 0.1666666666666667 Z[ATLAM5]
   + 0.1666666666666667 Z[LOSAng] + 0.1666666666666667 Z[DNVRng]
   + 0.1666666666666667 Z[WASHng] + 0.1666666666666667 Z[HSTNng]
Subject To
 flow-constr-t0-vKSCYng: X[IPLSng,KSCYng] + X[DNVRng,KSCYng]
   + X[HSTNng,KSCYng] - Y[KSCYng] = 0
 flow-constr-t0-vSNVAng: X[STTLng,SNVAng] + X[DNVRng,SNVAng]
   + X[LOSAng,SNVAng] - Y[SNVAng] = 0
 flow-constr-t0-vSTTLng: X[SNVAng,STTLng] + X[DNVRng,STTLng] - Y[STTLng]
   = 0
 flow-constr-t0-vATLAng: X[IPLSng,ATLAng] + X[WASHng,ATLAng]
   + X[ATLAM5,ATLAng] + X[HSTNng,ATLAng] - Y[ATLAng] = 0
 flow-constr-t0-vCHINng: X[IPLSng,CHINng] + X[NYCMng,CHINng] - Y[CHINng]
   = 0
 flow-constr-t0-vIPLSng: X[KSCYng,IPLSng] + X[CHINng,IPLSng]
   + X[ATLAng,IPLSng] - Y[IPLSng] = 0
 flow-constr-t0-vNYCMng: X[CHINng,NYCMng] + X[WASHng,NYCMng] - Y[NYCMng]
   = 0
 flow-constr-t0-vATLAM5: X[ATLAng,ATLAM5] - Y[ATLAM5] = 0
 flow-constr-t0-vLOSAng: X[SNVAng,LOSAng] + X[HSTNng,LOSAng] - Y[LOSAng]
   = 0
 flow-constr-t0-vDNVRng: X[SNVAng,DNVRng] + X[KSCYng,DNVRng]
   + X[STTLng,DNVRng] - Y[DNVRng] = 0
 flow-constr-t0-vWASHng: X[ATLAng,WASHng] + X[NYCMng,WASHng] - Y[WASHng]
   = 0
 flow-constr-t0-vHSTNng: X[LOSAng,HSTNng] + X[ATLAng,HSTNng]
   + X[KSCYng,HSTNng] - Y[HSTNng] = 0
 flow-constr-t1-vKSCYng: X[IPLSng,KSCYng] + X[DNVRng,KSCYng]
   + X[HSTNng,KSCYng] - Y[KSCYng] = 0
 flow-constr-t1-vSNVAng: X[STTLng,SNVAng] + X[DNVRng,SNVAng]
   + X[LOSAng,SNVAng] - Y[SNVAng] = 0
 flow-constr-t1-vSTTLng: X[SNVAng,STTLng] + X[DNVRng,STTLng] - Y[STTLng]
   = 0
 flow-constr-t1-vATLAng: X[IPLSng,ATLAng] + X[WASHng,ATLAng]
   + X[ATLAM5,ATLAng] + X[HSTNng,ATLAng] - Y[ATLAng] = 0
 flow-constr-t1-vCHINng: X[IPLSng,CHINng] + X[NYCMng,CHINng] - Y[CHINng]
   = 0
 flow-constr-t1-vIPLSng: X[KSCYng,IPLSng] + X[CHINng,IPLSng]
   + X[ATLAng,IPLSng] - Y[IPLSng] = 0
 flow-constr-t1-vNYCMng: X[CHINng,NYCMng] + X[WASHng,NYCMng] - Y[NYCMng]
   = 0
 flow-constr-t1-vATLAM5: X[ATLAng,ATLAM5] - Y[ATLAM5] = 0
 flow-constr-t1-vLOSAng: X[SNVAng,LOSAng] + X[HSTNng,LOSAng] - Y[LOSAng]
   = 0
 flow-constr-t1-vDNVRng: X[SNVAng,DNVRng] + X[KSCYng,DNVRng]
   + X[STTLng,DNVRng] - Y[DNVRng] = 0
 flow-constr-t1-vWASHng: X[ATLAng,WASHng] + X[NYCMng,WASHng] - Y[WASHng]
   = 0
 flow-constr-t1-vHSTNng: X[LOSAng,HSTNng] + X[ATLAng,HSTNng]
   + X[KSCYng,HSTNng] - Y[HSTNng] = 0
 flow-constr-t2-vKSCYng: X[IPLSng,KSCYng] + X[DNVRng,KSCYng]
   + X[HSTNng,KSCYng] - Y[KSCYng] = 0
 flow-constr-t2-vSNVAng: X[STTLng,SNVAng] + X[DNVRng,SNVAng]
   + X[LOSAng,SNVAng] - Y[SNVAng] = 0
 flow-constr-t2-vSTTLng: X[SNVAng,STTLng] + X[DNVRng,STTLng] - Y[STTLng]
   = 0
 flow-constr-t2-vATLAng: X[IPLSng,ATLAng] + X[WASHng,ATLAng]
   + X[ATLAM5,ATLAng] + X[HSTNng,ATLAng] - Y[ATLAng] = 0
 flow-constr-t2-vCHINng: X[IPLSng,CHINng] + X[NYCMng,CHINng] - Y[CHINng]
   = 0
 flow-constr-t2-vIPLSng: X[KSCYng,IPLSng] + X[CHINng,IPLSng]
   + X[ATLAng,IPLSng] - Y[IPLSng] = 0
 flow-constr-t2-vNYCMng: X[CHINng,NYCMng] + X[WASHng,NYCMng] - Y[NYCMng]
   = 0
 flow-constr-t2-vATLAM5: X[ATLAng,ATLAM5] - Y[ATLAM5] = 0
 flow-constr-t2-vLOSAng: X[SNVAng,LOSAng] + X[HSTNng,LOSAng] - Y[LOSAng]
   = 0
 flow-constr-t2-vDNVRng: X[SNVAng,DNVRng] + X[KSCYng,DNVRng]
   + X[STTLng,DNVRng] - Y[DNVRng] = 0
 flow-constr-t2-vWASHng: X[ATLAng,WASHng] + X[NYCMng,WASHng] - Y[WASHng]
   = 0
 flow-constr-t2-vHSTNng: X[LOSAng,HSTNng] + X[ATLAng,HSTNng]
   + X[KSCYng,HSTNng] - Y[HSTNng] = 0
 flow-constr-t3-vKSCYng: X[IPLSng,KSCYng] + X[DNVRng,KSCYng]
   + X[HSTNng,KSCYng] - Y[KSCYng] = 0
 flow-constr-t3-vSNVAng: X[STTLng,SNVAng] + X[DNVRng,SNVAng]
   + X[LOSAng,SNVAng] - Y[SNVAng] = 0
 flow-constr-t3-vSTTLng: X[SNVAng,STTLng] + X[DNVRng,STTLng] - Y[STTLng]
   = 0
 flow-constr-t3-vATLAng: X[IPLSng,ATLAng] + X[WASHng,ATLAng]
   + X[ATLAM5,ATLAng] + X[HSTNng,ATLAng] - Y[ATLAng] = 0
 flow-constr-t3-vCHINng: X[IPLSng,CHINng] + X[NYCMng,CHINng] - Y[CHINng]
   = 0
 flow-constr-t3-vIPLSng: X[KSCYng,IPLSng] + X[CHINng,IPLSng]
   + X[ATLAng,IPLSng] - Y[IPLSng] = 0
 flow-constr-t3-vNYCMng: X[CHINng,NYCMng] + X[WASHng,NYCMng] - Y[NYCMng]
   = 0
 flow-constr-t3-vATLAM5: X[ATLAng,ATLAM5] - Y[ATLAM5] = 0
 flow-constr-t3-vLOSAng: X[SNVAng,LOSAng] + X[HSTNng,LOSAng] - Y[LOSAng]
   = 0
 flow-constr-t3-vDNVRng: X[SNVAng,DNVRng] + X[KSCYng,DNVRng]
   + X[STTLng,DNVRng] - Y[DNVRng] = 0
 flow-constr-t3-vWASHng: X[ATLAng,WASHng] + X[NYCMng,WASHng] - Y[WASHng]
   = 0
 flow-constr-t3-vHSTNng: X[LOSAng,HSTNng] + X[ATLAng,HSTNng]
   + X[KSCYng,HSTNng] - Y[HSTNng] = 0
 flow-constr-t4-vKSCYng: X[IPLSng,KSCYng] + X[DNVRng,KSCYng]
   + X[HSTNng,KSCYng] - Y[KSCYng] = 0
 flow-constr-t4-vSNVAng: X[STTLng,SNVAng] + X[DNVRng,SNVAng]
   + X[LOSAng,SNVAng] - Y[SNVAng] = 0
 flow-constr-t4-vSTTLng: X[SNVAng,STTLng] + X[DNVRng,STTLng] - Y[STTLng]
   = 0
 flow-constr-t4-vATLAng: X[IPLSng,ATLAng] + X[WASHng,ATLAng]
   + X[ATLAM5,ATLAng] + X[HSTNng,ATLAng] - Y[ATLAng] = 0
 flow-constr-t4-vCHINng: X[IPLSng,CHINng] + X[NYCMng,CHINng] - Y[CHINng]
   = 0
 flow-constr-t4-vIPLSng: X[KSCYng,IPLSng] + X[CHINng,IPLSng]
   + X[ATLAng,IPLSng] - Y[IPLSng] = 0
 flow-constr-t4-vNYCMng: X[CHINng,NYCMng] + X[WASHng,NYCMng] - Y[NYCMng]
   = 0
 flow-constr-t4-vATLAM5: X[ATLAng,ATLAM5] - Y[ATLAM5] = 0
 flow-constr-t4-vLOSAng: X[SNVAng,LOSAng] + X[HSTNng,LOSAng] - Y[LOSAng]
   = 0
 flow-constr-t4-vDNVRng: X[SNVAng,DNVRng] + X[KSCYng,DNVRng]
   + X[STTLng,DNVRng] - Y[DNVRng] = 0
 flow-constr-t4-vWASHng: X[ATLAng,WASHng] + X[NYCMng,WASHng] - Y[WASHng]
   = 0
 flow-constr-t4-vHSTNng: X[LOSAng,HSTNng] + X[ATLAng,HSTNng]
   + X[KSCYng,HSTNng] - Y[HSTNng] = 0
 flow-constr-t5-vKSCYng: X[IPLSng,KSCYng] + X[DNVRng,KSCYng]
   + X[HSTNng,KSCYng] - Y[KSCYng] = 0
 flow-constr-t5-vSNVAng: X[STTLng,SNVAng] + X[DNVRng,SNVAng]
   + X[LOSAng,SNVAng] - Y[SNVAng] = 0
 flow-constr-t5-vSTTLng: X[SNVAng,STTLng] + X[DNVRng,STTLng] - Y[STTLng]
   = 0
 flow-constr-t5-vATLAng: X[IPLSng,ATLAng] + X[WASHng,ATLAng]
   + X[ATLAM5,ATLAng] + X[HSTNng,ATLAng] - Y[ATLAng] = 0
 flow-constr-t5-vCHINng: X[IPLSng,CHINng] + X[NYCMng,CHINng] - Y[CHINng]
   = 0
 flow-constr-t5-vIPLSng: X[KSCYng,IPLSng] + X[CHINng,IPLSng]
   + X[ATLAng,IPLSng] - Y[IPLSng] = 0
 flow-constr-t5-vNYCMng: X[CHINng,NYCMng] + X[WASHng,NYCMng] - Y[NYCMng]
   = 0
 flow-constr-t5-vATLAM5: X[ATLAng,ATLAM5] - Y[ATLAM5] = 0
 flow-constr-t5-vLOSAng: X[SNVAng,LOSAng] + X[HSTNng,LOSAng] - Y[LOSAng]
   = 0
 flow-constr-t5-vDNVRng: X[SNVAng,DNVRng] + X[KSCYng,DNVRng]
   + X[STTLng,DNVRng] - Y[DNVRng] = 0
 flow-constr-t5-vWASHng: X[ATLAng,WASHng] + X[NYCMng,WASHng] - Y[WASHng]
   = 0
 flow-constr-t5-vHSTNng: X[LOSAng,HSTNng] + X[ATLAng,HSTNng]
   + X[KSCYng,HSTNng] - Y[HSTNng] = 0
 root-constr-t0: Y[STTLng] + Y[LOSAng] + Y[DNVRng] >= 1
 root-constr-t1: Y[STTLng] + Y[LOSAng] + Y[DNVRng] >= 1
 root-constr-t2: Y[STTLng] + Y[LOSAng] + Y[DNVRng] >= 1
 root-constr-t3: Y[STTLng] + Y[LOSAng] + Y[DNVRng] >= 1
 root-constr-t4: Y[STTLng] + Y[LOSAng] + Y[DNVRng] >= 1
 root-constr-t5: Y[STTLng] + Y[LOSAng] + Y[DNVRng] >= 1
 weight-constr-t5: Z[KSCYng] + Z[KSCYng] + Z[KSCYng] + Z[KSCYng]
   + Z[KSCYng] + Z[KSCYng] >= 6
 weight-constr-t5: Z[STTLng] + Z[STTLng] + Z[STTLng] + Z[STTLng]
   + Z[STTLng] + Z[STTLng] >= 6
 weight-constr-t5: Z[ATLAng] + Z[ATLAng] + Z[ATLAng] + Z[ATLAng]
   + Z[ATLAng] + Z[ATLAng] >= 6
 weight-constr-t5: Z[CHINng] + Z[CHINng] + Z[CHINng] + Z[CHINng]
   + Z[CHINng] + Z[CHINng] >= 6
 weight-constr-t5: Z[IPLSng] + Z[IPLSng] + Z[IPLSng] + Z[IPLSng]
   + Z[IPLSng] + Z[IPLSng] >= 6
 weight-constr-t5: Z[NYCMng] + Z[NYCMng] + Z[NYCMng] + Z[NYCMng]
   + Z[NYCMng] + Z[NYCMng] >= 6
 weight-constr-t5: Z[ATLAM5] + Z[ATLAM5] + Z[ATLAM5] + Z[ATLAM5]
   + Z[ATLAM5] + Z[ATLAM5] >= 6
 weight-constr-t5: Z[LOSAng] + Z[LOSAng] + Z[LOSAng] + Z[LOSAng]
   + Z[LOSAng] + Z[LOSAng] >= 6
 weight-constr-t5: Z[DNVRng] + Z[DNVRng] + Z[DNVRng] + Z[DNVRng]
   + Z[DNVRng] + Z[DNVRng] >= 6
 weight-constr-t5: Z[WASHng] + Z[WASHng] + Z[WASHng] + Z[WASHng]
   + Z[WASHng] + Z[WASHng] >= 6
 weight-constr-t5: Z[HSTNng] + Z[HSTNng] + Z[HSTNng] + Z[HSTNng]
   + Z[HSTNng] + Z[HSTNng] >= 6
 strat-constr-t0-aCHINngNYCMng: - 1.048576e+06 X[CHINng,NYCMng] + Z[CHINng]
   - 2 Z[NYCMng] >= -1.048576e+06
 strat-constr-t0-aIPLSngATLAng: - 1.048576e+06 X[IPLSng,ATLAng]
   - 2 Z[ATLAng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t0-aIPLSngKSCYng: - 1.048576e+06 X[IPLSng,KSCYng]
   - 2 Z[KSCYng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t0-aKSCYngIPLSng: - 1.048576e+06 X[KSCYng,IPLSng] + Z[KSCYng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t0-aLOSAngHSTNng: - 1.048576e+06 X[LOSAng,HSTNng] + Z[LOSAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t0-aWASHngATLAng: - 1.048576e+06 X[WASHng,ATLAng]
   - 2 Z[ATLAng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t0-aCHINngIPLSng: - 1.048576e+06 X[CHINng,IPLSng] + Z[CHINng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t0-aIPLSngCHINng: - 1.048576e+06 X[IPLSng,CHINng]
   - 2 Z[CHINng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t0-aDNVRngSTTLng: - 1.048576e+06 X[DNVRng,STTLng]
   - 2 Z[STTLng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t0-aWASHngNYCMng: - 1.048576e+06 X[WASHng,NYCMng]
   - 2 Z[NYCMng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t0-aATLAM5ATLAng: - 1.048576e+06 X[ATLAM5,ATLAng]
   - 2 Z[ATLAng] + Z[ATLAM5] >= -1.048576e+06
 strat-constr-t0-aATLAngATLAM5: - 1.048576e+06 X[ATLAng,ATLAM5] + Z[ATLAng]
   - 2 Z[ATLAM5] >= -1.048576e+06
 strat-constr-t0-aATLAngIPLSng: - 1.048576e+06 X[ATLAng,IPLSng] + Z[ATLAng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t0-aDNVRngKSCYng: - 1.048576e+06 X[DNVRng,KSCYng]
   - 2 Z[KSCYng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t0-aATLAngHSTNng: - 1.048576e+06 X[ATLAng,HSTNng] + Z[ATLAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t0-aATLAngWASHng: - 1.048576e+06 X[ATLAng,WASHng] + Z[ATLAng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t0-aHSTNngKSCYng: - 1.048576e+06 X[HSTNng,KSCYng]
   - 2 Z[KSCYng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t0-aNYCMngWASHng: - 1.048576e+06 X[NYCMng,WASHng] + Z[NYCMng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t0-aHSTNngLOSAng: - 1.048576e+06 X[HSTNng,LOSAng]
   - 2 Z[LOSAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t0-aKSCYngDNVRng: - 1.048576e+06 X[KSCYng,DNVRng] + Z[KSCYng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t0-aNYCMngCHINng: - 1.048576e+06 X[NYCMng,CHINng]
   - 2 Z[CHINng] + Z[NYCMng] >= -1.048576e+06
 strat-constr-t0-aSTTLngDNVRng: - 1.048576e+06 X[STTLng,DNVRng] + Z[STTLng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t0-aHSTNngATLAng: - 1.048576e+06 X[HSTNng,ATLAng]
   - 2 Z[ATLAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t0-aKSCYngHSTNng: - 1.048576e+06 X[KSCYng,HSTNng] + Z[KSCYng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t1-aCHINngNYCMng: - 1.048576e+06 X[CHINng,NYCMng] + Z[CHINng]
   - 2 Z[NYCMng] >= -1.048576e+06
 strat-constr-t1-aIPLSngATLAng: - 1.048576e+06 X[IPLSng,ATLAng]
   - 2 Z[ATLAng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t1-aIPLSngKSCYng: - 1.048576e+06 X[IPLSng,KSCYng]
   - 2 Z[KSCYng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t1-aKSCYngIPLSng: - 1.048576e+06 X[KSCYng,IPLSng] + Z[KSCYng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t1-aLOSAngHSTNng: - 1.048576e+06 X[LOSAng,HSTNng] + Z[LOSAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t1-aWASHngATLAng: - 1.048576e+06 X[WASHng,ATLAng]
   - 2 Z[ATLAng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t1-aCHINngIPLSng: - 1.048576e+06 X[CHINng,IPLSng] + Z[CHINng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t1-aIPLSngCHINng: - 1.048576e+06 X[IPLSng,CHINng]
   - 2 Z[CHINng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t1-aDNVRngSTTLng: - 1.048576e+06 X[DNVRng,STTLng]
   - 2 Z[STTLng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t1-aWASHngNYCMng: - 1.048576e+06 X[WASHng,NYCMng]
   - 2 Z[NYCMng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t1-aATLAM5ATLAng: - 1.048576e+06 X[ATLAM5,ATLAng]
   - 2 Z[ATLAng] + Z[ATLAM5] >= -1.048576e+06
 strat-constr-t1-aATLAngATLAM5: - 1.048576e+06 X[ATLAng,ATLAM5] + Z[ATLAng]
   - 2 Z[ATLAM5] >= -1.048576e+06
 strat-constr-t1-aATLAngIPLSng: - 1.048576e+06 X[ATLAng,IPLSng] + Z[ATLAng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t1-aDNVRngKSCYng: - 1.048576e+06 X[DNVRng,KSCYng]
   - 2 Z[KSCYng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t1-aATLAngHSTNng: - 1.048576e+06 X[ATLAng,HSTNng] + Z[ATLAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t1-aATLAngWASHng: - 1.048576e+06 X[ATLAng,WASHng] + Z[ATLAng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t1-aHSTNngKSCYng: - 1.048576e+06 X[HSTNng,KSCYng]
   - 2 Z[KSCYng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t1-aNYCMngWASHng: - 1.048576e+06 X[NYCMng,WASHng] + Z[NYCMng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t1-aHSTNngLOSAng: - 1.048576e+06 X[HSTNng,LOSAng]
   - 2 Z[LOSAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t1-aKSCYngDNVRng: - 1.048576e+06 X[KSCYng,DNVRng] + Z[KSCYng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t1-aNYCMngCHINng: - 1.048576e+06 X[NYCMng,CHINng]
   - 2 Z[CHINng] + Z[NYCMng] >= -1.048576e+06
 strat-constr-t1-aSTTLngDNVRng: - 1.048576e+06 X[STTLng,DNVRng] + Z[STTLng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t1-aHSTNngATLAng: - 1.048576e+06 X[HSTNng,ATLAng]
   - 2 Z[ATLAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t1-aKSCYngHSTNng: - 1.048576e+06 X[KSCYng,HSTNng] + Z[KSCYng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t2-aCHINngNYCMng: - 1.048576e+06 X[CHINng,NYCMng] + Z[CHINng]
   - 2 Z[NYCMng] >= -1.048576e+06
 strat-constr-t2-aIPLSngATLAng: - 1.048576e+06 X[IPLSng,ATLAng]
   - 2 Z[ATLAng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t2-aIPLSngKSCYng: - 1.048576e+06 X[IPLSng,KSCYng]
   - 2 Z[KSCYng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t2-aKSCYngIPLSng: - 1.048576e+06 X[KSCYng,IPLSng] + Z[KSCYng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t2-aLOSAngHSTNng: - 1.048576e+06 X[LOSAng,HSTNng] + Z[LOSAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t2-aWASHngATLAng: - 1.048576e+06 X[WASHng,ATLAng]
   - 2 Z[ATLAng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t2-aCHINngIPLSng: - 1.048576e+06 X[CHINng,IPLSng] + Z[CHINng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t2-aIPLSngCHINng: - 1.048576e+06 X[IPLSng,CHINng]
   - 2 Z[CHINng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t2-aDNVRngSTTLng: - 1.048576e+06 X[DNVRng,STTLng]
   - 2 Z[STTLng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t2-aWASHngNYCMng: - 1.048576e+06 X[WASHng,NYCMng]
   - 2 Z[NYCMng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t2-aATLAM5ATLAng: - 1.048576e+06 X[ATLAM5,ATLAng]
   - 2 Z[ATLAng] + Z[ATLAM5] >= -1.048576e+06
 strat-constr-t2-aATLAngATLAM5: - 1.048576e+06 X[ATLAng,ATLAM5] + Z[ATLAng]
   - 2 Z[ATLAM5] >= -1.048576e+06
 strat-constr-t2-aATLAngIPLSng: - 1.048576e+06 X[ATLAng,IPLSng] + Z[ATLAng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t2-aDNVRngKSCYng: - 1.048576e+06 X[DNVRng,KSCYng]
   - 2 Z[KSCYng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t2-aATLAngHSTNng: - 1.048576e+06 X[ATLAng,HSTNng] + Z[ATLAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t2-aATLAngWASHng: - 1.048576e+06 X[ATLAng,WASHng] + Z[ATLAng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t2-aHSTNngKSCYng: - 1.048576e+06 X[HSTNng,KSCYng]
   - 2 Z[KSCYng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t2-aNYCMngWASHng: - 1.048576e+06 X[NYCMng,WASHng] + Z[NYCMng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t2-aHSTNngLOSAng: - 1.048576e+06 X[HSTNng,LOSAng]
   - 2 Z[LOSAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t2-aKSCYngDNVRng: - 1.048576e+06 X[KSCYng,DNVRng] + Z[KSCYng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t2-aNYCMngCHINng: - 1.048576e+06 X[NYCMng,CHINng]
   - 2 Z[CHINng] + Z[NYCMng] >= -1.048576e+06
 strat-constr-t2-aSTTLngDNVRng: - 1.048576e+06 X[STTLng,DNVRng] + Z[STTLng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t2-aHSTNngATLAng: - 1.048576e+06 X[HSTNng,ATLAng]
   - 2 Z[ATLAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t2-aKSCYngHSTNng: - 1.048576e+06 X[KSCYng,HSTNng] + Z[KSCYng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t3-aCHINngNYCMng: - 1.048576e+06 X[CHINng,NYCMng] + Z[CHINng]
   - 2 Z[NYCMng] >= -1.048576e+06
 strat-constr-t3-aIPLSngATLAng: - 1.048576e+06 X[IPLSng,ATLAng]
   - 2 Z[ATLAng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t3-aIPLSngKSCYng: - 1.048576e+06 X[IPLSng,KSCYng]
   - 2 Z[KSCYng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t3-aKSCYngIPLSng: - 1.048576e+06 X[KSCYng,IPLSng] + Z[KSCYng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t3-aLOSAngHSTNng: - 1.048576e+06 X[LOSAng,HSTNng] + Z[LOSAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t3-aWASHngATLAng: - 1.048576e+06 X[WASHng,ATLAng]
   - 2 Z[ATLAng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t3-aCHINngIPLSng: - 1.048576e+06 X[CHINng,IPLSng] + Z[CHINng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t3-aIPLSngCHINng: - 1.048576e+06 X[IPLSng,CHINng]
   - 2 Z[CHINng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t3-aDNVRngSTTLng: - 1.048576e+06 X[DNVRng,STTLng]
   - 2 Z[STTLng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t3-aWASHngNYCMng: - 1.048576e+06 X[WASHng,NYCMng]
   - 2 Z[NYCMng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t3-aATLAM5ATLAng: - 1.048576e+06 X[ATLAM5,ATLAng]
   - 2 Z[ATLAng] + Z[ATLAM5] >= -1.048576e+06
 strat-constr-t3-aATLAngATLAM5: - 1.048576e+06 X[ATLAng,ATLAM5] + Z[ATLAng]
   - 2 Z[ATLAM5] >= -1.048576e+06
 strat-constr-t3-aATLAngIPLSng: - 1.048576e+06 X[ATLAng,IPLSng] + Z[ATLAng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t3-aDNVRngKSCYng: - 1.048576e+06 X[DNVRng,KSCYng]
   - 2 Z[KSCYng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t3-aATLAngHSTNng: - 1.048576e+06 X[ATLAng,HSTNng] + Z[ATLAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t3-aATLAngWASHng: - 1.048576e+06 X[ATLAng,WASHng] + Z[ATLAng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t3-aHSTNngKSCYng: - 1.048576e+06 X[HSTNng,KSCYng]
   - 2 Z[KSCYng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t3-aNYCMngWASHng: - 1.048576e+06 X[NYCMng,WASHng] + Z[NYCMng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t3-aHSTNngLOSAng: - 1.048576e+06 X[HSTNng,LOSAng]
   - 2 Z[LOSAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t3-aKSCYngDNVRng: - 1.048576e+06 X[KSCYng,DNVRng] + Z[KSCYng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t3-aNYCMngCHINng: - 1.048576e+06 X[NYCMng,CHINng]
   - 2 Z[CHINng] + Z[NYCMng] >= -1.048576e+06
 strat-constr-t3-aSTTLngDNVRng: - 1.048576e+06 X[STTLng,DNVRng] + Z[STTLng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t3-aHSTNngATLAng: - 1.048576e+06 X[HSTNng,ATLAng]
   - 2 Z[ATLAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t3-aKSCYngHSTNng: - 1.048576e+06 X[KSCYng,HSTNng] + Z[KSCYng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t4-aCHINngNYCMng: - 1.048576e+06 X[CHINng,NYCMng] + Z[CHINng]
   - 2 Z[NYCMng] >= -1.048576e+06
 strat-constr-t4-aIPLSngATLAng: - 1.048576e+06 X[IPLSng,ATLAng]
   - 2 Z[ATLAng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t4-aIPLSngKSCYng: - 1.048576e+06 X[IPLSng,KSCYng]
   - 2 Z[KSCYng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t4-aKSCYngIPLSng: - 1.048576e+06 X[KSCYng,IPLSng] + Z[KSCYng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t4-aLOSAngHSTNng: - 1.048576e+06 X[LOSAng,HSTNng] + Z[LOSAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t4-aWASHngATLAng: - 1.048576e+06 X[WASHng,ATLAng]
   - 2 Z[ATLAng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t4-aCHINngIPLSng: - 1.048576e+06 X[CHINng,IPLSng] + Z[CHINng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t4-aIPLSngCHINng: - 1.048576e+06 X[IPLSng,CHINng]
   - 2 Z[CHINng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t4-aDNVRngSTTLng: - 1.048576e+06 X[DNVRng,STTLng]
   - 2 Z[STTLng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t4-aWASHngNYCMng: - 1.048576e+06 X[WASHng,NYCMng]
   - 2 Z[NYCMng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t4-aATLAM5ATLAng: - 1.048576e+06 X[ATLAM5,ATLAng]
   - 2 Z[ATLAng] + Z[ATLAM5] >= -1.048576e+06
 strat-constr-t4-aATLAngATLAM5: - 1.048576e+06 X[ATLAng,ATLAM5] + Z[ATLAng]
   - 2 Z[ATLAM5] >= -1.048576e+06
 strat-constr-t4-aATLAngIPLSng: - 1.048576e+06 X[ATLAng,IPLSng] + Z[ATLAng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t4-aDNVRngKSCYng: - 1.048576e+06 X[DNVRng,KSCYng]
   - 2 Z[KSCYng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t4-aATLAngHSTNng: - 1.048576e+06 X[ATLAng,HSTNng] + Z[ATLAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t4-aATLAngWASHng: - 1.048576e+06 X[ATLAng,WASHng] + Z[ATLAng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t4-aHSTNngKSCYng: - 1.048576e+06 X[HSTNng,KSCYng]
   - 2 Z[KSCYng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t4-aNYCMngWASHng: - 1.048576e+06 X[NYCMng,WASHng] + Z[NYCMng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t4-aHSTNngLOSAng: - 1.048576e+06 X[HSTNng,LOSAng]
   - 2 Z[LOSAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t4-aKSCYngDNVRng: - 1.048576e+06 X[KSCYng,DNVRng] + Z[KSCYng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t4-aNYCMngCHINng: - 1.048576e+06 X[NYCMng,CHINng]
   - 2 Z[CHINng] + Z[NYCMng] >= -1.048576e+06
 strat-constr-t4-aSTTLngDNVRng: - 1.048576e+06 X[STTLng,DNVRng] + Z[STTLng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t4-aHSTNngATLAng: - 1.048576e+06 X[HSTNng,ATLAng]
   - 2 Z[ATLAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t4-aKSCYngHSTNng: - 1.048576e+06 X[KSCYng,HSTNng] + Z[KSCYng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t5-aCHINngNYCMng: - 1.048576e+06 X[CHINng,NYCMng] + Z[CHINng]
   - 2 Z[NYCMng] >= -1.048576e+06
 strat-constr-t5-aIPLSngATLAng: - 1.048576e+06 X[IPLSng,ATLAng]
   - 2 Z[ATLAng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t5-aIPLSngKSCYng: - 1.048576e+06 X[IPLSng,KSCYng]
   - 2 Z[KSCYng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t5-aKSCYngIPLSng: - 1.048576e+06 X[KSCYng,IPLSng] + Z[KSCYng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t5-aLOSAngHSTNng: - 1.048576e+06 X[LOSAng,HSTNng] + Z[LOSAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t5-aWASHngATLAng: - 1.048576e+06 X[WASHng,ATLAng]
   - 2 Z[ATLAng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t5-aCHINngIPLSng: - 1.048576e+06 X[CHINng,IPLSng] + Z[CHINng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t5-aIPLSngCHINng: - 1.048576e+06 X[IPLSng,CHINng]
   - 2 Z[CHINng] + Z[IPLSng] >= -1.048576e+06
 strat-constr-t5-aDNVRngSTTLng: - 1.048576e+06 X[DNVRng,STTLng]
   - 2 Z[STTLng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t5-aWASHngNYCMng: - 1.048576e+06 X[WASHng,NYCMng]
   - 2 Z[NYCMng] + Z[WASHng] >= -1.048576e+06
 strat-constr-t5-aATLAM5ATLAng: - 1.048576e+06 X[ATLAM5,ATLAng]
   - 2 Z[ATLAng] + Z[ATLAM5] >= -1.048576e+06
 strat-constr-t5-aATLAngATLAM5: - 1.048576e+06 X[ATLAng,ATLAM5] + Z[ATLAng]
   - 2 Z[ATLAM5] >= -1.048576e+06
 strat-constr-t5-aATLAngIPLSng: - 1.048576e+06 X[ATLAng,IPLSng] + Z[ATLAng]
   - 2 Z[IPLSng] >= -1.048576e+06
 strat-constr-t5-aDNVRngKSCYng: - 1.048576e+06 X[DNVRng,KSCYng]
   - 2 Z[KSCYng] + Z[DNVRng] >= -1.048576e+06
 strat-constr-t5-aATLAngHSTNng: - 1.048576e+06 X[ATLAng,HSTNng] + Z[ATLAng]
   - 2 Z[HSTNng] >= -1.048576e+06
 strat-constr-t5-aATLAngWASHng: - 1.048576e+06 X[ATLAng,WASHng] + Z[ATLAng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t5-aHSTNngKSCYng: - 1.048576e+06 X[HSTNng,KSCYng]
   - 2 Z[KSCYng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t5-aNYCMngWASHng: - 1.048576e+06 X[NYCMng,WASHng] + Z[NYCMng]
   - 2 Z[WASHng] >= -1.048576e+06
 strat-constr-t5-aHSTNngLOSAng: - 1.048576e+06 X[HSTNng,LOSAng]
   - 2 Z[LOSAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t5-aKSCYngDNVRng: - 1.048576e+06 X[KSCYng,DNVRng] + Z[KSCYng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t5-aNYCMngCHINng: - 1.048576e+06 X[NYCMng,CHINng]
   - 2 Z[CHINng] + Z[NYCMng] >= -1.048576e+06
 strat-constr-t5-aSTTLngDNVRng: - 1.048576e+06 X[STTLng,DNVRng] + Z[STTLng]
   - 2 Z[DNVRng] >= -1.048576e+06
 strat-constr-t5-aHSTNngATLAng: - 1.048576e+06 X[HSTNng,ATLAng]
   - 2 Z[ATLAng] + Z[HSTNng] >= -1.048576e+06
 strat-constr-t5-aKSCYngHSTNng: - 1.048576e+06 X[KSCYng,HSTNng] + Z[KSCYng]
   - 2 Z[HSTNng] >= -1.048576e+06
 weight-constr2-t0-vKSCYng: - 524288 Y[KSCYng] + Z[KSCYng] <= 0
 weight-constr3-t0-vKSCYng: Y[KSCYng] - Z[KSCYng] <= 0
 weight-constr2-t0-vSNVAng: - 524288 Y[SNVAng] + Z[SNVAng] <= 0
 weight-constr3-t0-vSNVAng: Y[SNVAng] - Z[SNVAng] <= 0
 weight-constr2-t0-vSTTLng: - 524288 Y[STTLng] + Z[STTLng] <= 0
 weight-constr3-t0-vSTTLng: Y[STTLng] - Z[STTLng] <= 0
 weight-constr2-t0-vATLAng: - 524288 Y[ATLAng] + Z[ATLAng] <= 0
 weight-constr3-t0-vATLAng: Y[ATLAng] - Z[ATLAng] <= 0
 weight-constr2-t0-vCHINng: - 524288 Y[CHINng] + Z[CHINng] <= 0
 weight-constr3-t0-vCHINng: Y[CHINng] - Z[CHINng] <= 0
 weight-constr2-t0-vIPLSng: - 524288 Y[IPLSng] + Z[IPLSng] <= 0
 weight-constr3-t0-vIPLSng: Y[IPLSng] - Z[IPLSng] <= 0
 weight-constr2-t0-vNYCMng: - 524288 Y[NYCMng] + Z[NYCMng] <= 0
 weight-constr3-t0-vNYCMng: Y[NYCMng] - Z[NYCMng] <= 0
 weight-constr2-t0-vATLAM5: - 524288 Y[ATLAM5] + Z[ATLAM5] <= 0
 weight-constr3-t0-vATLAM5: Y[ATLAM5] - Z[ATLAM5] <= 0
 weight-constr2-t0-vLOSAng: - 524288 Y[LOSAng] + Z[LOSAng] <= 0
 weight-constr3-t0-vLOSAng: Y[LOSAng] - Z[LOSAng] <= 0
 weight-constr2-t0-vDNVRng: - 524288 Y[DNVRng] + Z[DNVRng] <= 0
 weight-constr3-t0-vDNVRng: Y[DNVRng] - Z[DNVRng] <= 0
 weight-constr2-t0-vWASHng: - 524288 Y[WASHng] + Z[WASHng] <= 0
 weight-constr3-t0-vWASHng: Y[WASHng] - Z[WASHng] <= 0
 weight-constr2-t0-vHSTNng: - 524288 Y[HSTNng] + Z[HSTNng] <= 0
 weight-constr3-t0-vHSTNng: Y[HSTNng] - Z[HSTNng] <= 0
 weight-constr2-t1-vKSCYng: - 524288 Y[KSCYng] + Z[KSCYng] <= 0
 weight-constr3-t1-vKSCYng: Y[KSCYng] - Z[KSCYng] <= 0
 weight-constr2-t1-vSNVAng: - 524288 Y[SNVAng] + Z[SNVAng] <= 0
 weight-constr3-t1-vSNVAng: Y[SNVAng] - Z[SNVAng] <= 0
 weight-constr2-t1-vSTTLng: - 524288 Y[STTLng] + Z[STTLng] <= 0
 weight-constr3-t1-vSTTLng: Y[STTLng] - Z[STTLng] <= 0
 weight-constr2-t1-vATLAng: - 524288 Y[ATLAng] + Z[ATLAng] <= 0
 weight-constr3-t1-vATLAng: Y[ATLAng] - Z[ATLAng] <= 0
 weight-constr2-t1-vCHINng: - 524288 Y[CHINng] + Z[CHINng] <= 0
 weight-constr3-t1-vCHINng: Y[CHINng] - Z[CHINng] <= 0
 weight-constr2-t1-vIPLSng: - 524288 Y[IPLSng] + Z[IPLSng] <= 0
 weight-constr3-t1-vIPLSng: Y[IPLSng] - Z[IPLSng] <= 0
 weight-constr2-t1-vNYCMng: - 524288 Y[NYCMng] + Z[NYCMng] <= 0
 weight-constr3-t1-vNYCMng: Y[NYCMng] - Z[NYCMng] <= 0
 weight-constr2-t1-vATLAM5: - 524288 Y[ATLAM5] + Z[ATLAM5] <= 0
 weight-constr3-t1-vATLAM5: Y[ATLAM5] - Z[ATLAM5] <= 0
 weight-constr2-t1-vLOSAng: - 524288 Y[LOSAng] + Z[LOSAng] <= 0
 weight-constr3-t1-vLOSAng: Y[LOSAng] - Z[LOSAng] <= 0
 weight-constr2-t1-vDNVRng: - 524288 Y[DNVRng] + Z[DNVRng] <= 0
 weight-constr3-t1-vDNVRng: Y[DNVRng] - Z[DNVRng] <= 0
 weight-constr2-t1-vWASHng: - 524288 Y[WASHng] + Z[WASHng] <= 0
 weight-constr3-t1-vWASHng: Y[WASHng] - Z[WASHng] <= 0
 weight-constr2-t1-vHSTNng: - 524288 Y[HSTNng] + Z[HSTNng] <= 0
 weight-constr3-t1-vHSTNng: Y[HSTNng] - Z[HSTNng] <= 0
 weight-constr2-t2-vKSCYng: - 524288 Y[KSCYng] + Z[KSCYng] <= 0
 weight-constr3-t2-vKSCYng: Y[KSCYng] - Z[KSCYng] <= 0
 weight-constr2-t2-vSNVAng: - 524288 Y[SNVAng] + Z[SNVAng] <= 0
 weight-constr3-t2-vSNVAng: Y[SNVAng] - Z[SNVAng] <= 0
 weight-constr2-t2-vSTTLng: - 524288 Y[STTLng] + Z[STTLng] <= 0
 weight-constr3-t2-vSTTLng: Y[STTLng] - Z[STTLng] <= 0
 weight-constr2-t2-vATLAng: - 524288 Y[ATLAng] + Z[ATLAng] <= 0
 weight-constr3-t2-vATLAng: Y[ATLAng] - Z[ATLAng] <= 0
 weight-constr2-t2-vCHINng: - 524288 Y[CHINng] + Z[CHINng] <= 0
 weight-constr3-t2-vCHINng: Y[CHINng] - Z[CHINng] <= 0
 weight-constr2-t2-vIPLSng: - 524288 Y[IPLSng] + Z[IPLSng] <= 0
 weight-constr3-t2-vIPLSng: Y[IPLSng] - Z[IPLSng] <= 0
 weight-constr2-t2-vNYCMng: - 524288 Y[NYCMng] + Z[NYCMng] <= 0
 weight-constr3-t2-vNYCMng: Y[NYCMng] - Z[NYCMng] <= 0
 weight-constr2-t2-vATLAM5: - 524288 Y[ATLAM5] + Z[ATLAM5] <= 0
 weight-constr3-t2-vATLAM5: Y[ATLAM5] - Z[ATLAM5] <= 0
 weight-constr2-t2-vLOSAng: - 524288 Y[LOSAng] + Z[LOSAng] <= 0
 weight-constr3-t2-vLOSAng: Y[LOSAng] - Z[LOSAng] <= 0
 weight-constr2-t2-vDNVRng: - 524288 Y[DNVRng] + Z[DNVRng] <= 0
 weight-constr3-t2-vDNVRng: Y[DNVRng] - Z[DNVRng] <= 0
 weight-constr2-t2-vWASHng: - 524288 Y[WASHng] + Z[WASHng] <= 0
 weight-constr3-t2-vWASHng: Y[WASHng] - Z[WASHng] <= 0
 weight-constr2-t2-vHSTNng: - 524288 Y[HSTNng] + Z[HSTNng] <= 0
 weight-constr3-t2-vHSTNng: Y[HSTNng] - Z[HSTNng] <= 0
 weight-constr2-t3-vKSCYng: - 524288 Y[KSCYng] + Z[KSCYng] <= 0
 weight-constr3-t3-vKSCYng: Y[KSCYng] - Z[KSCYng] <= 0
 weight-constr2-t3-vSNVAng: - 524288 Y[SNVAng] + Z[SNVAng] <= 0
 weight-constr3-t3-vSNVAng: Y[SNVAng] - Z[SNVAng] <= 0
 weight-constr2-t3-vSTTLng: - 524288 Y[STTLng] + Z[STTLng] <= 0
 weight-constr3-t3-vSTTLng: Y[STTLng] - Z[STTLng] <= 0
 weight-constr2-t3-vATLAng: - 524288 Y[ATLAng] + Z[ATLAng] <= 0
 weight-constr3-t3-vATLAng: Y[ATLAng] - Z[ATLAng] <= 0
 weight-constr2-t3-vCHINng: - 524288 Y[CHINng] + Z[CHINng] <= 0
 weight-constr3-t3-vCHINng: Y[CHINng] - Z[CHINng] <= 0
 weight-constr2-t3-vIPLSng: - 524288 Y[IPLSng] + Z[IPLSng] <= 0
 weight-constr3-t3-vIPLSng: Y[IPLSng] - Z[IPLSng] <= 0
 weight-constr2-t3-vNYCMng: - 524288 Y[NYCMng] + Z[NYCMng] <= 0
 weight-constr3-t3-vNYCMng: Y[NYCMng] - Z[NYCMng] <= 0
 weight-constr2-t3-vATLAM5: - 524288 Y[ATLAM5] + Z[ATLAM5] <= 0
 weight-constr3-t3-vATLAM5: Y[ATLAM5] - Z[ATLAM5] <= 0
 weight-constr2-t3-vLOSAng: - 524288 Y[LOSAng] + Z[LOSAng] <= 0
 weight-constr3-t3-vLOSAng: Y[LOSAng] - Z[LOSAng] <= 0
 weight-constr2-t3-vDNVRng: - 524288 Y[DNVRng] + Z[DNVRng] <= 0
 weight-constr3-t3-vDNVRng: Y[DNVRng] - Z[DNVRng] <= 0
 weight-constr2-t3-vWASHng: - 524288 Y[WASHng] + Z[WASHng] <= 0
 weight-constr3-t3-vWASHng: Y[WASHng] - Z[WASHng] <= 0
 weight-constr2-t3-vHSTNng: - 524288 Y[HSTNng] + Z[HSTNng] <= 0
 weight-constr3-t3-vHSTNng: Y[HSTNng] - Z[HSTNng] <= 0
 weight-constr2-t4-vKSCYng: - 524288 Y[KSCYng] + Z[KSCYng] <= 0
 weight-constr3-t4-vKSCYng: Y[KSCYng] - Z[KSCYng] <= 0
 weight-constr2-t4-vSNVAng: - 524288 Y[SNVAng] + Z[SNVAng] <= 0
 weight-constr3-t4-vSNVAng: Y[SNVAng] - Z[SNVAng] <= 0
 weight-constr2-t4-vSTTLng: - 524288 Y[STTLng] + Z[STTLng] <= 0
 weight-constr3-t4-vSTTLng: Y[STTLng] - Z[STTLng] <= 0
 weight-constr2-t4-vATLAng: - 524288 Y[ATLAng] + Z[ATLAng] <= 0
 weight-constr3-t4-vATLAng: Y[ATLAng] - Z[ATLAng] <= 0
 weight-constr2-t4-vCHINng: - 524288 Y[CHINng] + Z[CHINng] <= 0
 weight-constr3-t4-vCHINng: Y[CHINng] - Z[CHINng] <= 0
 weight-constr2-t4-vIPLSng: - 524288 Y[IPLSng] + Z[IPLSng] <= 0
 weight-constr3-t4-vIPLSng: Y[IPLSng] - Z[IPLSng] <= 0
 weight-constr2-t4-vNYCMng: - 524288 Y[NYCMng] + Z[NYCMng] <= 0
 weight-constr3-t4-vNYCMng: Y[NYCMng] - Z[NYCMng] <= 0
 weight-constr2-t4-vATLAM5: - 524288 Y[ATLAM5] + Z[ATLAM5] <= 0
 weight-constr3-t4-vATLAM5: Y[ATLAM5] - Z[ATLAM5] <= 0
 weight-constr2-t4-vLOSAng: - 524288 Y[LOSAng] + Z[LOSAng] <= 0
 weight-constr3-t4-vLOSAng: Y[LOSAng] - Z[LOSAng] <= 0
 weight-constr2-t4-vDNVRng: - 524288 Y[DNVRng] + Z[DNVRng] <= 0
 weight-constr3-t4-vDNVRng: Y[DNVRng] - Z[DNVRng] <= 0
 weight-constr2-t4-vWASHng: - 524288 Y[WASHng] + Z[WASHng] <= 0
 weight-constr3-t4-vWASHng: Y[WASHng] - Z[WASHng] <= 0
 weight-constr2-t4-vHSTNng: - 524288 Y[HSTNng] + Z[HSTNng] <= 0
 weight-constr3-t4-vHSTNng: Y[HSTNng] - Z[HSTNng] <= 0
 weight-constr2-t5-vKSCYng: - 524288 Y[KSCYng] + Z[KSCYng] <= 0
 weight-constr3-t5-vKSCYng: Y[KSCYng] - Z[KSCYng] <= 0
 weight-constr2-t5-vSNVAng: - 524288 Y[SNVAng] + Z[SNVAng] <= 0
 weight-constr3-t5-vSNVAng: Y[SNVAng] - Z[SNVAng] <= 0
 weight-constr2-t5-vSTTLng: - 524288 Y[STTLng] + Z[STTLng] <= 0
 weight-constr3-t5-vSTTLng: Y[STTLng] - Z[STTLng] <= 0
 weight-constr2-t5-vATLAng: - 524288 Y[ATLAng] + Z[ATLAng] <= 0
 weight-constr3-t5-vATLAng: Y[ATLAng] - Z[ATLAng] <= 0
 weight-constr2-t5-vCHINng: - 524288 Y[CHINng] + Z[CHINng] <= 0
 weight-constr3-t5-vCHINng: Y[CHINng] - Z[CHINng] <= 0
 weight-constr2-t5-vIPLSng: - 524288 Y[IPLSng] + Z[IPLSng] <= 0
 weight-constr3-t5-vIPLSng: Y[IPLSng] - Z[IPLSng] <= 0
 weight-constr2-t5-vNYCMng: - 524288 Y[NYCMng] + Z[NYCMng] <= 0
 weight-constr3-t5-vNYCMng: Y[NYCMng] - Z[NYCMng] <= 0
 weight-constr2-t5-vATLAM5: - 524288 Y[ATLAM5] + Z[ATLAM5] <= 0
 weight-constr3-t5-vATLAM5: Y[ATLAM5] - Z[ATLAM5] <= 0
 weight-constr2-t5-vLOSAng: - 524288 Y[LOSAng] + Z[LOSAng] <= 0
 weight-constr3-t5-vLOSAng: Y[LOSAng] - Z[LOSAng] <= 0
 weight-constr2-t5-vDNVRng: - 524288 Y[DNVRng] + Z[DNVRng] <= 0
 weight-constr3-t5-vDNVRng: Y[DNVRng] - Z[DNVRng] <= 0
 weight-constr2-t5-vWASHng: - 524288 Y[WASHng] + Z[WASHng] <= 0
 weight-constr3-t5-vWASHng: Y[WASHng] - Z[WASHng] <= 0
 weight-constr2-t5-vHSTNng: - 524288 Y[HSTNng] + Z[HSTNng] <= 0
 weight-constr3-t5-vHSTNng: Y[HSTNng] - Z[HSTNng] <= 0
Bounds
 Z[KSCYng] <= 524288
 Z[SNVAng] <= 524288
 Z[STTLng] <= 524288
 Z[ATLAng] <= 524288
 Z[CHINng] <= 524288
 Z[IPLSng] <= 524288
 Z[NYCMng] <= 524288
 Z[ATLAM5] <= 524288
 Z[LOSAng] <= 524288
 Z[DNVRng] <= 524288
 Z[WASHng] <= 524288
 Z[HSTNng] <= 524288
 Z[KSCYng] <= 524288
 Z[SNVAng] <= 524288
 Z[STTLng] <= 524288
 Z[ATLAng] <= 524288
 Z[CHINng] <= 524288
 Z[IPLSng] <= 524288
 Z[NYCMng] <= 524288
 Z[ATLAM5] <= 524288
 Z[LOSAng] <= 524288
 Z[DNVRng] <= 524288
 Z[WASHng] <= 524288
 Z[HSTNng] <= 524288
 Z[KSCYng] <= 524288
 Z[SNVAng] <= 524288
 Z[STTLng] <= 524288
 Z[ATLAng] <= 524288
 Z[CHINng] <= 524288
 Z[IPLSng] <= 524288
 Z[NYCMng] <= 524288
 Z[ATLAM5] <= 524288
 Z[LOSAng] <= 524288
 Z[DNVRng] <= 524288
 Z[WASHng] <= 524288
 Z[HSTNng] <= 524288
 Z[KSCYng] <= 524288
 Z[SNVAng] <= 524288
 Z[STTLng] <= 524288
 Z[ATLAng] <= 524288
 Z[CHINng] <= 524288
 Z[IPLSng] <= 524288
 Z[NYCMng] <= 524288
 Z[ATLAM5] <= 524288
 Z[LOSAng] <= 524288
 Z[DNVRng] <= 524288
 Z[WASHng] <= 524288
 Z[HSTNng] <= 524288
 Z[KSCYng] <= 524288
 Z[SNVAng] <= 524288
 Z[STTLng] <= 524288
 Z[ATLAng] <= 524288
 Z[CHINng] <= 524288
 Z[IPLSng] <= 524288
 Z[NYCMng] <= 524288
 Z[ATLAM5] <= 524288
 Z[LOSAng] <= 524288
 Z[DNVRng] <= 524288
 Z[WASHng] <= 524288
 Z[HSTNng] <= 524288
 Z[KSCYng] <= 524288
 Z[SNVAng] <= 524288
 Z[STTLng] <= 524288
 Z[ATLAng] <= 524288
 Z[CHINng] <= 524288
 Z[IPLSng] <= 524288
 Z[NYCMng] <= 524288
 Z[ATLAM5] <= 524288
 Z[LOSAng] <= 524288
 Z[DNVRng] <= 524288
 Z[WASHng] <= 524288
 Z[HSTNng] <= 524288
Binaries
 X[SNVAng,STTLng] X[CHINng,NYCMng] X[IPLSng,ATLAng] X[IPLSng,KSCYng]
 X[KSCYng,IPLSng] X[LOSAng,HSTNng] X[WASHng,ATLAng] X[CHINng,IPLSng]
 X[IPLSng,CHINng] X[DNVRng,STTLng] X[WASHng,NYCMng] X[ATLAM5,ATLAng]
 X[STTLng,SNVAng] X[ATLAng,ATLAM5] X[ATLAng,IPLSng] X[SNVAng,LOSAng]
 X[DNVRng,KSCYng] X[ATLAng,HSTNng] X[DNVRng,SNVAng] X[SNVAng,DNVRng]
 X[ATLAng,WASHng] X[HSTNng,KSCYng] X[NYCMng,WASHng] X[HSTNng,LOSAng]
 X[LOSAng,SNVAng] X[KSCYng,DNVRng] X[NYCMng,CHINng] X[STTLng,DNVRng]
 X[HSTNng,ATLAng] X[KSCYng,HSTNng] Y[KSCYng] Y[SNVAng] Y[STTLng] Y[ATLAng]
 Y[CHINng] Y[IPLSng] Y[NYCMng] Y[ATLAM5] Y[LOSAng] Y[DNVRng] Y[WASHng]
 Y[HSTNng] X[SNVAng,STTLng] X[CHINng,NYCMng] X[IPLSng,ATLAng]
 X[IPLSng,KSCYng] X[KSCYng,IPLSng] X[LOSAng,HSTNng] X[WASHng,ATLAng]
 X[CHINng,IPLSng] X[IPLSng,CHINng] X[DNVRng,STTLng] X[WASHng,NYCMng]
 X[ATLAM5,ATLAng] X[STTLng,SNVAng] X[ATLAng,ATLAM5] X[ATLAng,IPLSng]
 X[SNVAng,LOSAng] X[DNVRng,KSCYng] X[ATLAng,HSTNng] X[DNVRng,SNVAng]
 X[SNVAng,DNVRng] X[ATLAng,WASHng] X[HSTNng,KSCYng] X[NYCMng,WASHng]
 X[HSTNng,LOSAng] X[LOSAng,SNVAng] X[KSCYng,DNVRng] X[NYCMng,CHINng]
 X[STTLng,DNVRng] X[HSTNng,ATLAng] X[KSCYng,HSTNng] Y[KSCYng] Y[SNVAng]
 Y[STTLng] Y[ATLAng] Y[CHINng] Y[IPLSng] Y[NYCMng] Y[ATLAM5] Y[LOSAng]
 Y[DNVRng] Y[WASHng] Y[HSTNng] X[SNVAng,STTLng] X[CHINng,NYCMng]
 X[IPLSng,ATLAng] X[IPLSng,KSCYng] X[KSCYng,IPLSng] X[LOSAng,HSTNng]
 X[WASHng,ATLAng] X[CHINng,IPLSng] X[IPLSng,CHINng] X[DNVRng,STTLng]
 X[WASHng,NYCMng] X[ATLAM5,ATLAng] X[STTLng,SNVAng] X[ATLAng,ATLAM5]
 X[ATLAng,IPLSng] X[SNVAng,LOSAng] X[DNVRng,KSCYng] X[ATLAng,HSTNng]
 X[DNVRng,SNVAng] X[SNVAng,DNVRng] X[ATLAng,WASHng] X[HSTNng,KSCYng]
 X[NYCMng,WASHng] X[HSTNng,LOSAng] X[LOSAng,SNVAng] X[KSCYng,DNVRng]
 X[NYCMng,CHINng] X[STTLng,DNVRng] X[HSTNng,ATLAng] X[KSCYng,HSTNng]
 Y[KSCYng] Y[SNVAng] Y[STTLng] Y[ATLAng] Y[CHINng] Y[IPLSng] Y[NYCMng]
 Y[ATLAM5] Y[LOSAng] Y[DNVRng] Y[WASHng] Y[HSTNng] X[SNVAng,STTLng]
 X[CHINng,NYCMng] X[IPLSng,ATLAng] X[IPLSng,KSCYng] X[KSCYng,IPLSng]
 X[LOSAng,HSTNng] X[WASHng,ATLAng] X[CHINng,IPLSng] X[IPLSng,CHINng]
 X[DNVRng,STTLng] X[WASHng,NYCMng] X[ATLAM5,ATLAng] X[STTLng,SNVAng]
 X[ATLAng,ATLAM5] X[ATLAng,IPLSng] X[SNVAng,LOSAng] X[DNVRng,KSCYng]
 X[ATLAng,HSTNng] X[DNVRng,SNVAng] X[SNVAng,DNVRng] X[ATLAng,WASHng]
 X[HSTNng,KSCYng] X[NYCMng,WASHng] X[HSTNng,LOSAng] X[LOSAng,SNVAng]
 X[KSCYng,DNVRng] X[NYCMng,CHINng] X[STTLng,DNVRng] X[HSTNng,ATLAng]
 X[KSCYng,HSTNng] Y[KSCYng] Y[SNVAng] Y[STTLng] Y[ATLAng] Y[CHINng]
 Y[IPLSng] Y[NYCMng] Y[ATLAM5] Y[LOSAng] Y[DNVRng] Y[WASHng] Y[HSTNng]
 X[SNVAng,STTLng] X[CHINng,NYCMng] X[IPLSng,ATLAng] X[IPLSng,KSCYng]
 X[KSCYng,IPLSng] X[LOSAng,HSTNng] X[WASHng,ATLAng] X[CHINng,IPLSng]
 X[IPLSng,CHINng] X[DNVRng,STTLng] X[WASHng,NYCMng] X[ATLAM5,ATLAng]
 X[STTLng,SNVAng] X[ATLAng,ATLAM5] X[ATLAng,IPLSng] X[SNVAng,LOSAng]
 X[DNVRng,KSCYng] X[ATLAng,HSTNng] X[DNVRng,SNVAng] X[SNVAng,DNVRng]
 X[ATLAng,WASHng] X[HSTNng,KSCYng] X[NYCMng,WASHng] X[HSTNng,LOSAng]
 X[LOSAng,SNVAng] X[KSCYng,DNVRng] X[NYCMng,CHINng] X[STTLng,DNVRng]
 X[HSTNng,ATLAng] X[KSCYng,HSTNng] Y[KSCYng] Y[SNVAng] Y[STTLng] Y[ATLAng]
 Y[CHINng] Y[IPLSng] Y[NYCMng] Y[ATLAM5] Y[LOSAng] Y[DNVRng] Y[WASHng]
 Y[HSTNng] X[SNVAng,STTLng] X[CHINng,NYCMng] X[IPLSng,ATLAng]
 X[IPLSng,KSCYng] X[KSCYng,IPLSng] X[LOSAng,HSTNng] X[WASHng,ATLAng]
 X[CHINng,IPLSng] X[IPLSng,CHINng] X[DNVRng,STTLng] X[WASHng,NYCMng]
 X[ATLAM5,ATLAng] X[STTLng,SNVAng] X[ATLAng,ATLAM5] X[ATLAng,IPLSng]
 X[SNVAng,LOSAng] X[DNVRng,KSCYng] X[ATLAng,HSTNng] X[DNVRng,SNVAng]
 X[SNVAng,DNVRng] X[ATLAng,WASHng] X[HSTNng,KSCYng] X[NYCMng,WASHng]
 X[HSTNng,LOSAng] X[LOSAng,SNVAng] X[KSCYng,DNVRng] X[NYCMng,CHINng]
 X[STTLng,DNVRng] X[HSTNng,ATLAng] X[KSCYng,HSTNng] Y[KSCYng] Y[SNVAng]
 Y[STTLng] Y[ATLAng] Y[CHINng] Y[IPLSng] Y[NYCMng] Y[ATLAM5] Y[LOSAng]
 Y[DNVRng] Y[WASHng] Y[HSTNng]
End
