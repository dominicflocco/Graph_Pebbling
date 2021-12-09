\ Model weight-function-optimizer
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  0 V[1] + V[0] + V[3] + V[2] + V[5] + V[4] + V[7] + V[6]
Subject To
 R0: 12 V[0] + 6 V[2] + 6 V[4] + 3 V[6] <= 27
 R1: 12 V[3] + 6 V[2] + 6 V[7] + 3 V[6] <= 27
 R2: 12 V[5] + 6 V[4] + 6 V[7] + 3 V[6] <= 27
Bounds
 V[1] <= 256
 V[0] <= 256
 V[3] <= 256
 V[2] <= 256
 V[5] <= 256
 V[4] <= 256
 V[7] <= 256
 V[6] <= 256
Generals
 V[1] V[0] V[3] V[2] V[5] V[4] V[7] V[6]
End
