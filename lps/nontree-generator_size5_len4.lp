\ Model nontree-generator
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0.2 Z[1] + 0.2 Z[0] + 0.2 Z[3] + 0.2 Z[2] + 0.2 Z[5] + 0.2 Z[4] + 0.2 Z[7]
   + 0.2 Z[6] + 0.2 Z[1] + 0.2 Z[0] + 0.2 Z[3] + 0.2 Z[2] + 0.2 Z[5]
   + 0.2 Z[4] + 0.2 Z[7] + 0.2 Z[6] + 0.2 Z[1] + 0.2 Z[0] + 0.2 Z[3]
   + 0.2 Z[2] + 0.2 Z[5] + 0.2 Z[4] + 0.2 Z[7] + 0.2 Z[6] + 0.2 Z[1]
   + 0.2 Z[0] + 0.2 Z[3] + 0.2 Z[2] + 0.2 Z[5] + 0.2 Z[4] + 0.2 Z[7]
   + 0.2 Z[6] + 0.2 Z[1] + 0.2 Z[0] + 0.2 Z[3] + 0.2 Z[2] + 0.2 Z[5]
   + 0.2 Z[4] + 0.2 Z[7] + 0.2 Z[6]
Subject To
 flow-constr-t0-v1: X[0,1] + X[3,1] - Y[1] = 0
 flow-constr-t0-v0: X[1,0] + X[2,0] - Y[0] = 0
 flow-constr-t0-v3: X[1,3] + X[4,3] + X[5,3] + X[6,3] + X[7,3] - Y[3] = 0
 flow-constr-t0-v2: X[0,2] + X[4,2] + X[5,2] + X[6,2] - Y[2] = 0
 flow-constr-t0-v5: X[2,5] + X[3,5] + X[7,5] - Y[5] = 0
 flow-constr-t0-v4: X[2,4] + X[3,4] + X[7,4] - Y[4] = 0
 flow-constr-t0-v7: X[3,7] + X[4,7] + X[5,7] + X[6,7] - Y[7] = 0
 flow-constr-t0-v6: X[2,6] + X[3,6] + X[7,6] - Y[6] = 0
 flow-constr-t1-v1: X[0,1] + X[3,1] - Y[1] = 0
 flow-constr-t1-v0: X[1,0] + X[2,0] - Y[0] = 0
 flow-constr-t1-v3: X[1,3] + X[4,3] + X[5,3] + X[6,3] + X[7,3] - Y[3] = 0
 flow-constr-t1-v2: X[0,2] + X[4,2] + X[5,2] + X[6,2] - Y[2] = 0
 flow-constr-t1-v5: X[2,5] + X[3,5] + X[7,5] - Y[5] = 0
 flow-constr-t1-v4: X[2,4] + X[3,4] + X[7,4] - Y[4] = 0
 flow-constr-t1-v7: X[3,7] + X[4,7] + X[5,7] + X[6,7] - Y[7] = 0
 flow-constr-t1-v6: X[2,6] + X[3,6] + X[7,6] - Y[6] = 0
 flow-constr-t2-v1: X[0,1] + X[3,1] - Y[1] = 0
 flow-constr-t2-v0: X[1,0] + X[2,0] - Y[0] = 0
 flow-constr-t2-v3: X[1,3] + X[4,3] + X[5,3] + X[6,3] + X[7,3] - Y[3] = 0
 flow-constr-t2-v2: X[0,2] + X[4,2] + X[5,2] + X[6,2] - Y[2] = 0
 flow-constr-t2-v5: X[2,5] + X[3,5] + X[7,5] - Y[5] = 0
 flow-constr-t2-v4: X[2,4] + X[3,4] + X[7,4] - Y[4] = 0
 flow-constr-t2-v7: X[3,7] + X[4,7] + X[5,7] + X[6,7] - Y[7] = 0
 flow-constr-t2-v6: X[2,6] + X[3,6] + X[7,6] - Y[6] = 0
 flow-constr-t3-v1: X[0,1] + X[3,1] - Y[1] = 0
 flow-constr-t3-v0: X[1,0] + X[2,0] - Y[0] = 0
 flow-constr-t3-v3: X[1,3] + X[4,3] + X[5,3] + X[6,3] + X[7,3] - Y[3] = 0
 flow-constr-t3-v2: X[0,2] + X[4,2] + X[5,2] + X[6,2] - Y[2] = 0
 flow-constr-t3-v5: X[2,5] + X[3,5] + X[7,5] - Y[5] = 0
 flow-constr-t3-v4: X[2,4] + X[3,4] + X[7,4] - Y[4] = 0
 flow-constr-t3-v7: X[3,7] + X[4,7] + X[5,7] + X[6,7] - Y[7] = 0
 flow-constr-t3-v6: X[2,6] + X[3,6] + X[7,6] - Y[6] = 0
 flow-constr-t4-v1: X[0,1] + X[3,1] - Y[1] = 0
 flow-constr-t4-v0: X[1,0] + X[2,0] - Y[0] = 0
 flow-constr-t4-v3: X[1,3] + X[4,3] + X[5,3] + X[6,3] + X[7,3] - Y[3] = 0
 flow-constr-t4-v2: X[0,2] + X[4,2] + X[5,2] + X[6,2] - Y[2] = 0
 flow-constr-t4-v5: X[2,5] + X[3,5] + X[7,5] - Y[5] = 0
 flow-constr-t4-v4: X[2,4] + X[3,4] + X[7,4] - Y[4] = 0
 flow-constr-t4-v7: X[3,7] + X[4,7] + X[5,7] + X[6,7] - Y[7] = 0
 flow-constr-t4-v6: X[2,6] + X[3,6] + X[7,6] - Y[6] = 0
 root-constr-t0: Y[3] + Y[2] + Y[7] = 1
 root-constr-t1: Y[3] + Y[2] + Y[7] = 1
 root-constr-t2: Y[3] + Y[2] + Y[7] = 1
 root-constr-t3: Y[3] + Y[2] + Y[7] = 1
 root-constr-t4: Y[3] + Y[2] + Y[7] = 1
 weight-constr-t4: Z[1] + Z[1] + Z[1] + Z[1] + Z[1] >= 5
 weight-constr-t4: Z[0] + Z[0] + Z[0] + Z[0] + Z[0] >= 5
 weight-constr-t4: Z[3] + Z[3] + Z[3] + Z[3] + Z[3] >= 5
 weight-constr-t4: Z[2] + Z[2] + Z[2] + Z[2] + Z[2] >= 5
 weight-constr-t4: Z[5] + Z[5] + Z[5] + Z[5] + Z[5] >= 5
 weight-constr-t4: Z[4] + Z[4] + Z[4] + Z[4] + Z[4] >= 5
 weight-constr-t4: Z[7] + Z[7] + Z[7] + Z[7] + Z[7] >= 5
 strat-constr-t0-a10: - 16 X[1,0] + Z[1] - 2 Z[0] >= -16
 strat-constr-t0-a01: - 16 X[0,1] - 2 Z[1] + Z[0] >= -16
 strat-constr-t0-a20: - 16 X[2,0] - 2 Z[0] + Z[2] >= -16
 strat-constr-t0-a02: - 16 X[0,2] + Z[0] - 2 Z[2] >= -16
 strat-constr-t0-a31: - 16 X[3,1] - 2 Z[1] + Z[3] >= -16
 strat-constr-t0-a13: - 16 X[1,3] + Z[1] - 2 Z[3] >= -16
 strat-constr-t0-a42: - 16 X[4,2] - 2 Z[2] + Z[4] >= -16
 strat-constr-t0-a24: - 16 X[2,4] + Z[2] - 2 Z[4] >= -16
 strat-constr-t0-a52: - 16 X[5,2] - 2 Z[2] + Z[5] >= -16
 strat-constr-t0-a25: - 16 X[2,5] + Z[2] - 2 Z[5] >= -16
 strat-constr-t0-a43: - 16 X[4,3] - 2 Z[3] + Z[4] >= -16
 strat-constr-t0-a34: - 16 X[3,4] + Z[3] - 2 Z[4] >= -16
 strat-constr-t0-a53: - 16 X[5,3] - 2 Z[3] + Z[5] >= -16
 strat-constr-t0-a35: - 16 X[3,5] + Z[3] - 2 Z[5] >= -16
 strat-constr-t0-a73: - 16 X[7,3] - 2 Z[3] + Z[7] >= -16
 strat-constr-t0-a37: - 16 X[3,7] + Z[3] - 2 Z[7] >= -16
 strat-constr-t0-a74: - 16 X[7,4] - 2 Z[4] + Z[7] >= -16
 strat-constr-t0-a47: - 16 X[4,7] + Z[4] - 2 Z[7] >= -16
 strat-constr-t0-a75: - 16 X[7,5] - 2 Z[5] + Z[7] >= -16
 strat-constr-t0-a57: - 16 X[5,7] + Z[5] - 2 Z[7] >= -16
 strat-constr-t1-a10: - 16 X[1,0] + Z[1] - 2 Z[0] >= -16
 strat-constr-t1-a01: - 16 X[0,1] - 2 Z[1] + Z[0] >= -16
 strat-constr-t1-a20: - 16 X[2,0] - 2 Z[0] + Z[2] >= -16
 strat-constr-t1-a02: - 16 X[0,2] + Z[0] - 2 Z[2] >= -16
 strat-constr-t1-a31: - 16 X[3,1] - 2 Z[1] + Z[3] >= -16
 strat-constr-t1-a13: - 16 X[1,3] + Z[1] - 2 Z[3] >= -16
 strat-constr-t1-a42: - 16 X[4,2] - 2 Z[2] + Z[4] >= -16
 strat-constr-t1-a24: - 16 X[2,4] + Z[2] - 2 Z[4] >= -16
 strat-constr-t1-a52: - 16 X[5,2] - 2 Z[2] + Z[5] >= -16
 strat-constr-t1-a25: - 16 X[2,5] + Z[2] - 2 Z[5] >= -16
 strat-constr-t1-a43: - 16 X[4,3] - 2 Z[3] + Z[4] >= -16
 strat-constr-t1-a34: - 16 X[3,4] + Z[3] - 2 Z[4] >= -16
 strat-constr-t1-a53: - 16 X[5,3] - 2 Z[3] + Z[5] >= -16
 strat-constr-t1-a35: - 16 X[3,5] + Z[3] - 2 Z[5] >= -16
 strat-constr-t1-a73: - 16 X[7,3] - 2 Z[3] + Z[7] >= -16
 strat-constr-t1-a37: - 16 X[3,7] + Z[3] - 2 Z[7] >= -16
 strat-constr-t1-a74: - 16 X[7,4] - 2 Z[4] + Z[7] >= -16
 strat-constr-t1-a47: - 16 X[4,7] + Z[4] - 2 Z[7] >= -16
 strat-constr-t1-a75: - 16 X[7,5] - 2 Z[5] + Z[7] >= -16
 strat-constr-t1-a57: - 16 X[5,7] + Z[5] - 2 Z[7] >= -16
 strat-constr-t2-a10: - 16 X[1,0] + Z[1] - 2 Z[0] >= -16
 strat-constr-t2-a01: - 16 X[0,1] - 2 Z[1] + Z[0] >= -16
 strat-constr-t2-a20: - 16 X[2,0] - 2 Z[0] + Z[2] >= -16
 strat-constr-t2-a02: - 16 X[0,2] + Z[0] - 2 Z[2] >= -16
 strat-constr-t2-a31: - 16 X[3,1] - 2 Z[1] + Z[3] >= -16
 strat-constr-t2-a13: - 16 X[1,3] + Z[1] - 2 Z[3] >= -16
 strat-constr-t2-a42: - 16 X[4,2] - 2 Z[2] + Z[4] >= -16
 strat-constr-t2-a24: - 16 X[2,4] + Z[2] - 2 Z[4] >= -16
 strat-constr-t2-a52: - 16 X[5,2] - 2 Z[2] + Z[5] >= -16
 strat-constr-t2-a25: - 16 X[2,5] + Z[2] - 2 Z[5] >= -16
 strat-constr-t2-a43: - 16 X[4,3] - 2 Z[3] + Z[4] >= -16
 strat-constr-t2-a34: - 16 X[3,4] + Z[3] - 2 Z[4] >= -16
 strat-constr-t2-a53: - 16 X[5,3] - 2 Z[3] + Z[5] >= -16
 strat-constr-t2-a35: - 16 X[3,5] + Z[3] - 2 Z[5] >= -16
 strat-constr-t2-a73: - 16 X[7,3] - 2 Z[3] + Z[7] >= -16
 strat-constr-t2-a37: - 16 X[3,7] + Z[3] - 2 Z[7] >= -16
 strat-constr-t2-a74: - 16 X[7,4] - 2 Z[4] + Z[7] >= -16
 strat-constr-t2-a47: - 16 X[4,7] + Z[4] - 2 Z[7] >= -16
 strat-constr-t2-a75: - 16 X[7,5] - 2 Z[5] + Z[7] >= -16
 strat-constr-t2-a57: - 16 X[5,7] + Z[5] - 2 Z[7] >= -16
 strat-constr-t3-a10: - 16 X[1,0] + Z[1] - 2 Z[0] >= -16
 strat-constr-t3-a01: - 16 X[0,1] - 2 Z[1] + Z[0] >= -16
 strat-constr-t3-a20: - 16 X[2,0] - 2 Z[0] + Z[2] >= -16
 strat-constr-t3-a02: - 16 X[0,2] + Z[0] - 2 Z[2] >= -16
 strat-constr-t3-a31: - 16 X[3,1] - 2 Z[1] + Z[3] >= -16
 strat-constr-t3-a13: - 16 X[1,3] + Z[1] - 2 Z[3] >= -16
 strat-constr-t3-a42: - 16 X[4,2] - 2 Z[2] + Z[4] >= -16
 strat-constr-t3-a24: - 16 X[2,4] + Z[2] - 2 Z[4] >= -16
 strat-constr-t3-a52: - 16 X[5,2] - 2 Z[2] + Z[5] >= -16
 strat-constr-t3-a25: - 16 X[2,5] + Z[2] - 2 Z[5] >= -16
 strat-constr-t3-a43: - 16 X[4,3] - 2 Z[3] + Z[4] >= -16
 strat-constr-t3-a34: - 16 X[3,4] + Z[3] - 2 Z[4] >= -16
 strat-constr-t3-a53: - 16 X[5,3] - 2 Z[3] + Z[5] >= -16
 strat-constr-t3-a35: - 16 X[3,5] + Z[3] - 2 Z[5] >= -16
 strat-constr-t3-a73: - 16 X[7,3] - 2 Z[3] + Z[7] >= -16
 strat-constr-t3-a37: - 16 X[3,7] + Z[3] - 2 Z[7] >= -16
 strat-constr-t3-a74: - 16 X[7,4] - 2 Z[4] + Z[7] >= -16
 strat-constr-t3-a47: - 16 X[4,7] + Z[4] - 2 Z[7] >= -16
 strat-constr-t3-a75: - 16 X[7,5] - 2 Z[5] + Z[7] >= -16
 strat-constr-t3-a57: - 16 X[5,7] + Z[5] - 2 Z[7] >= -16
 strat-constr-t4-a10: - 16 X[1,0] + Z[1] - 2 Z[0] >= -16
 strat-constr-t4-a01: - 16 X[0,1] - 2 Z[1] + Z[0] >= -16
 strat-constr-t4-a20: - 16 X[2,0] - 2 Z[0] + Z[2] >= -16
 strat-constr-t4-a02: - 16 X[0,2] + Z[0] - 2 Z[2] >= -16
 strat-constr-t4-a31: - 16 X[3,1] - 2 Z[1] + Z[3] >= -16
 strat-constr-t4-a13: - 16 X[1,3] + Z[1] - 2 Z[3] >= -16
 strat-constr-t4-a42: - 16 X[4,2] - 2 Z[2] + Z[4] >= -16
 strat-constr-t4-a24: - 16 X[2,4] + Z[2] - 2 Z[4] >= -16
 strat-constr-t4-a52: - 16 X[5,2] - 2 Z[2] + Z[5] >= -16
 strat-constr-t4-a25: - 16 X[2,5] + Z[2] - 2 Z[5] >= -16
 strat-constr-t4-a43: - 16 X[4,3] - 2 Z[3] + Z[4] >= -16
 strat-constr-t4-a34: - 16 X[3,4] + Z[3] - 2 Z[4] >= -16
 strat-constr-t4-a53: - 16 X[5,3] - 2 Z[3] + Z[5] >= -16
 strat-constr-t4-a35: - 16 X[3,5] + Z[3] - 2 Z[5] >= -16
 strat-constr-t4-a73: - 16 X[7,3] - 2 Z[3] + Z[7] >= -16
 strat-constr-t4-a37: - 16 X[3,7] + Z[3] - 2 Z[7] >= -16
 strat-constr-t4-a74: - 16 X[7,4] - 2 Z[4] + Z[7] >= -16
 strat-constr-t4-a47: - 16 X[4,7] + Z[4] - 2 Z[7] >= -16
 strat-constr-t4-a75: - 16 X[7,5] - 2 Z[5] + Z[7] >= -16
 strat-constr-t4-a57: - 16 X[5,7] + Z[5] - 2 Z[7] >= -16
 R152: 2 X[1,0] + 2 X[0,1] + 2 X[2,0] + 2 X[0,2] + 2 X[3,1] + 2 X[1,3]
   + 2 X[4,2] + 2 X[2,4] + 2 X[5,2] + 2 X[2,5] + 2 X[4,3] + 2 X[3,4]
   + 2 X[5,3] + 2 X[3,5] + 2 X[7,3] + 2 X[3,7] + 2 X[7,4] + 2 X[4,7]
   + 2 X[7,5] + 2 X[5,7] - 2 Y[1] - 2 Y[0] - 2 Y[3] - 2 Y[2] - 2 Y[5]
   - 2 Y[4] - 2 Y[7] - 2 Y[6] = -1
 R153: - X[6,2] - X[2,6] - X[6,3] - X[3,6] - X[7,6] - X[6,7] = -1
 R154: 2 X[1,0] + 2 X[0,1] + 2 X[2,0] + 2 X[0,2] + 2 X[3,1] + 2 X[1,3]
   + 2 X[4,2] + 2 X[2,4] + 2 X[5,2] + 2 X[2,5] + 2 X[4,3] + 2 X[3,4]
   + 2 X[5,3] + 2 X[3,5] + 2 X[7,3] + 2 X[3,7] + 2 X[7,4] + 2 X[4,7]
   + 2 X[7,5] + 2 X[5,7] - 2 Y[1] - 2 Y[0] - 2 Y[3] - 2 Y[2] - 2 Y[5]
   - 2 Y[4] - 2 Y[7] - 2 Y[6] = -1
 R155: - X[6,2] - X[2,6] - X[6,3] - X[3,6] - X[7,6] - X[6,7] = -1
 R156: 2 X[1,0] + 2 X[0,1] + 2 X[2,0] + 2 X[0,2] + 2 X[3,1] + 2 X[1,3]
   + 2 X[4,2] + 2 X[2,4] + 2 X[5,2] + 2 X[2,5] + 2 X[4,3] + 2 X[3,4]
   + 2 X[5,3] + 2 X[3,5] + 2 X[7,3] + 2 X[3,7] + 2 X[7,4] + 2 X[4,7]
   + 2 X[7,5] + 2 X[5,7] - 2 Y[1] - 2 Y[0] - 2 Y[3] - 2 Y[2] - 2 Y[5]
   - 2 Y[4] - 2 Y[7] - 2 Y[6] = -1
 R157: - X[6,2] - X[2,6] - X[6,3] - X[3,6] - X[7,6] - X[6,7] = -1
 R158: 2 X[1,0] + 2 X[0,1] + 2 X[2,0] + 2 X[0,2] + 2 X[3,1] + 2 X[1,3]
   + 2 X[4,2] + 2 X[2,4] + 2 X[5,2] + 2 X[2,5] + 2 X[4,3] + 2 X[3,4]
   + 2 X[5,3] + 2 X[3,5] + 2 X[7,3] + 2 X[3,7] + 2 X[7,4] + 2 X[4,7]
   + 2 X[7,5] + 2 X[5,7] - 2 Y[1] - 2 Y[0] - 2 Y[3] - 2 Y[2] - 2 Y[5]
   - 2 Y[4] - 2 Y[7] - 2 Y[6] = -1
 R159: - X[6,2] - X[2,6] - X[6,3] - X[3,6] - X[7,6] - X[6,7] = -1
 R160: 2 X[1,0] + 2 X[0,1] + 2 X[2,0] + 2 X[0,2] + 2 X[3,1] + 2 X[1,3]
   + 2 X[4,2] + 2 X[2,4] + 2 X[5,2] + 2 X[2,5] + 2 X[4,3] + 2 X[3,4]
   + 2 X[5,3] + 2 X[3,5] + 2 X[7,3] + 2 X[3,7] + 2 X[7,4] + 2 X[4,7]
   + 2 X[7,5] + 2 X[5,7] - 2 Y[1] - 2 Y[0] - 2 Y[3] - 2 Y[2] - 2 Y[5]
   - 2 Y[4] - 2 Y[7] - 2 Y[6] = -1
 R161: - X[6,2] - X[2,6] - X[6,3] - X[3,6] - X[7,6] - X[6,7] = -1
 weight-constr2-t0-v1: - 8 Y[1] + Z[1] <= 0
 weight-constr3-t0-v1: Y[1] - Z[1] <= 0
 weight-constr2-t0-v0: - 8 Y[0] + Z[0] <= 0
 weight-constr3-t0-v0: Y[0] - Z[0] <= 0
 weight-constr2-t0-v3: - 8 Y[3] + Z[3] <= 0
 weight-constr3-t0-v3: Y[3] - Z[3] <= 0
 weight-constr2-t0-v2: - 8 Y[2] + Z[2] <= 0
 weight-constr3-t0-v2: Y[2] - Z[2] <= 0
 weight-constr2-t0-v5: - 8 Y[5] + Z[5] <= 0
 weight-constr3-t0-v5: Y[5] - Z[5] <= 0
 weight-constr2-t0-v4: - 8 Y[4] + Z[4] <= 0
 weight-constr3-t0-v4: Y[4] - Z[4] <= 0
 weight-constr2-t0-v7: - 8 Y[7] + Z[7] <= 0
 weight-constr3-t0-v7: Y[7] - Z[7] <= 0
 weight-constr2-t0-v6: - 8 Y[6] + Z[6] <= 0
 weight-constr3-t0-v6: Y[6] - Z[6] <= 0
 weight-constr2-t1-v1: - 8 Y[1] + Z[1] <= 0
 weight-constr3-t1-v1: Y[1] - Z[1] <= 0
 weight-constr2-t1-v0: - 8 Y[0] + Z[0] <= 0
 weight-constr3-t1-v0: Y[0] - Z[0] <= 0
 weight-constr2-t1-v3: - 8 Y[3] + Z[3] <= 0
 weight-constr3-t1-v3: Y[3] - Z[3] <= 0
 weight-constr2-t1-v2: - 8 Y[2] + Z[2] <= 0
 weight-constr3-t1-v2: Y[2] - Z[2] <= 0
 weight-constr2-t1-v5: - 8 Y[5] + Z[5] <= 0
 weight-constr3-t1-v5: Y[5] - Z[5] <= 0
 weight-constr2-t1-v4: - 8 Y[4] + Z[4] <= 0
 weight-constr3-t1-v4: Y[4] - Z[4] <= 0
 weight-constr2-t1-v7: - 8 Y[7] + Z[7] <= 0
 weight-constr3-t1-v7: Y[7] - Z[7] <= 0
 weight-constr2-t1-v6: - 8 Y[6] + Z[6] <= 0
 weight-constr3-t1-v6: Y[6] - Z[6] <= 0
 weight-constr2-t2-v1: - 8 Y[1] + Z[1] <= 0
 weight-constr3-t2-v1: Y[1] - Z[1] <= 0
 weight-constr2-t2-v0: - 8 Y[0] + Z[0] <= 0
 weight-constr3-t2-v0: Y[0] - Z[0] <= 0
 weight-constr2-t2-v3: - 8 Y[3] + Z[3] <= 0
 weight-constr3-t2-v3: Y[3] - Z[3] <= 0
 weight-constr2-t2-v2: - 8 Y[2] + Z[2] <= 0
 weight-constr3-t2-v2: Y[2] - Z[2] <= 0
 weight-constr2-t2-v5: - 8 Y[5] + Z[5] <= 0
 weight-constr3-t2-v5: Y[5] - Z[5] <= 0
 weight-constr2-t2-v4: - 8 Y[4] + Z[4] <= 0
 weight-constr3-t2-v4: Y[4] - Z[4] <= 0
 weight-constr2-t2-v7: - 8 Y[7] + Z[7] <= 0
 weight-constr3-t2-v7: Y[7] - Z[7] <= 0
 weight-constr2-t2-v6: - 8 Y[6] + Z[6] <= 0
 weight-constr3-t2-v6: Y[6] - Z[6] <= 0
 weight-constr2-t3-v1: - 8 Y[1] + Z[1] <= 0
 weight-constr3-t3-v1: Y[1] - Z[1] <= 0
 weight-constr2-t3-v0: - 8 Y[0] + Z[0] <= 0
 weight-constr3-t3-v0: Y[0] - Z[0] <= 0
 weight-constr2-t3-v3: - 8 Y[3] + Z[3] <= 0
 weight-constr3-t3-v3: Y[3] - Z[3] <= 0
 weight-constr2-t3-v2: - 8 Y[2] + Z[2] <= 0
 weight-constr3-t3-v2: Y[2] - Z[2] <= 0
 weight-constr2-t3-v5: - 8 Y[5] + Z[5] <= 0
 weight-constr3-t3-v5: Y[5] - Z[5] <= 0
 weight-constr2-t3-v4: - 8 Y[4] + Z[4] <= 0
 weight-constr3-t3-v4: Y[4] - Z[4] <= 0
 weight-constr2-t3-v7: - 8 Y[7] + Z[7] <= 0
 weight-constr3-t3-v7: Y[7] - Z[7] <= 0
 weight-constr2-t3-v6: - 8 Y[6] + Z[6] <= 0
 weight-constr3-t3-v6: Y[6] - Z[6] <= 0
 weight-constr2-t4-v1: - 8 Y[1] + Z[1] <= 0
 weight-constr3-t4-v1: Y[1] - Z[1] <= 0
 weight-constr2-t4-v0: - 8 Y[0] + Z[0] <= 0
 weight-constr3-t4-v0: Y[0] - Z[0] <= 0
 weight-constr2-t4-v3: - 8 Y[3] + Z[3] <= 0
 weight-constr3-t4-v3: Y[3] - Z[3] <= 0
 weight-constr2-t4-v2: - 8 Y[2] + Z[2] <= 0
 weight-constr3-t4-v2: Y[2] - Z[2] <= 0
 weight-constr2-t4-v5: - 8 Y[5] + Z[5] <= 0
 weight-constr3-t4-v5: Y[5] - Z[5] <= 0
 weight-constr2-t4-v4: - 8 Y[4] + Z[4] <= 0
 weight-constr3-t4-v4: Y[4] - Z[4] <= 0
 weight-constr2-t4-v7: - 8 Y[7] + Z[7] <= 0
 weight-constr3-t4-v7: Y[7] - Z[7] <= 0
 weight-constr2-t4-v6: - 8 Y[6] + Z[6] <= 0
 weight-constr3-t4-v6: Y[6] - Z[6] <= 0
Bounds
 Z[1] <= 8
 Z[0] <= 8
 Z[3] <= 8
 Z[2] <= 8
 Z[5] <= 8
 Z[4] <= 8
 Z[7] <= 8
 Z[6] <= 8
 Z[1] <= 8
 Z[0] <= 8
 Z[3] <= 8
 Z[2] <= 8
 Z[5] <= 8
 Z[4] <= 8
 Z[7] <= 8
 Z[6] <= 8
 Z[1] <= 8
 Z[0] <= 8
 Z[3] <= 8
 Z[2] <= 8
 Z[5] <= 8
 Z[4] <= 8
 Z[7] <= 8
 Z[6] <= 8
 Z[1] <= 8
 Z[0] <= 8
 Z[3] <= 8
 Z[2] <= 8
 Z[5] <= 8
 Z[4] <= 8
 Z[7] <= 8
 Z[6] <= 8
 Z[1] <= 8
 Z[0] <= 8
 Z[3] <= 8
 Z[2] <= 8
 Z[5] <= 8
 Z[4] <= 8
 Z[7] <= 8
 Z[6] <= 8
Binaries
 X[1,0] X[0,1] X[2,0] X[0,2] X[3,1] X[1,3] X[4,2] X[2,4] X[5,2] X[2,5]
 X[6,2] X[2,6] X[4,3] X[3,4] X[5,3] X[3,5] X[6,3] X[3,6] X[7,3] X[3,7]
 X[7,4] X[4,7] X[7,5] X[5,7] X[7,6] X[6,7] Y[1] Y[0] Y[3] Y[2] Y[5] Y[4]
 Y[7] Y[6] X[1,0] X[0,1] X[2,0] X[0,2] X[3,1] X[1,3] X[4,2] X[2,4] X[5,2]
 X[2,5] X[6,2] X[2,6] X[4,3] X[3,4] X[5,3] X[3,5] X[6,3] X[3,6] X[7,3]
 X[3,7] X[7,4] X[4,7] X[7,5] X[5,7] X[7,6] X[6,7] Y[1] Y[0] Y[3] Y[2] Y[5]
 Y[4] Y[7] Y[6] X[1,0] X[0,1] X[2,0] X[0,2] X[3,1] X[1,3] X[4,2] X[2,4]
 X[5,2] X[2,5] X[6,2] X[2,6] X[4,3] X[3,4] X[5,3] X[3,5] X[6,3] X[3,6]
 X[7,3] X[3,7] X[7,4] X[4,7] X[7,5] X[5,7] X[7,6] X[6,7] Y[1] Y[0] Y[3]
 Y[2] Y[5] Y[4] Y[7] Y[6] X[1,0] X[0,1] X[2,0] X[0,2] X[3,1] X[1,3] X[4,2]
 X[2,4] X[5,2] X[2,5] X[6,2] X[2,6] X[4,3] X[3,4] X[5,3] X[3,5] X[6,3]
 X[3,6] X[7,3] X[3,7] X[7,4] X[4,7] X[7,5] X[5,7] X[7,6] X[6,7] Y[1] Y[0]
 Y[3] Y[2] Y[5] Y[4] Y[7] Y[6] X[1,0] X[0,1] X[2,0] X[0,2] X[3,1] X[1,3]
 X[4,2] X[2,4] X[5,2] X[2,5] X[6,2] X[2,6] X[4,3] X[3,4] X[5,3] X[3,5]
 X[6,3] X[3,6] X[7,3] X[3,7] X[7,4] X[4,7] X[7,5] X[5,7] X[7,6] X[6,7] Y[1]
 Y[0] Y[3] Y[2] Y[5] Y[4] Y[7] Y[6]
Generals
 Z[1] Z[0] Z[3] Z[2] Z[5] Z[4] Z[7] Z[6] Z[1] Z[0] Z[3] Z[2] Z[5] Z[4] Z[7]
 Z[6] Z[1] Z[0] Z[3] Z[2] Z[5] Z[4] Z[7] Z[6] Z[1] Z[0] Z[3] Z[2] Z[5] Z[4]
 Z[7] Z[6] Z[1] Z[0] Z[3] Z[2] Z[5] Z[4] Z[7] Z[6]
End
