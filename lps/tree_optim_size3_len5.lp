\ Model tree-optimizer-symmetry
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0.3333333333333333 Z[(0,_1)] + 0.3333333333333333 Z[(0,_0)]
   + 0.3333333333333333 Z[(1,_0)] + 0.3333333333333333 Z[(1,_1)]
   + 0.3333333333333333 Z[(2,_1)] + 0.3333333333333333 Z[(1,_2)]
   + 0.3333333333333333 Z[(2,_0)] + 0.3333333333333333 Z[(0,_2)]
   + 0.3333333333333333 Z[(2,_2)] + 0.3333333333333333 Z[(0,_1)]
   + 0.3333333333333333 Z[(0,_0)] + 0.3333333333333333 Z[(1,_0)]
   + 0.3333333333333333 Z[(1,_1)] + 0.3333333333333333 Z[(2,_1)]
   + 0.3333333333333333 Z[(1,_2)] + 0.3333333333333333 Z[(2,_0)]
   + 0.3333333333333333 Z[(0,_2)] + 0.3333333333333333 Z[(2,_2)]
   + 0.3333333333333333 Z[(0,_1)] + 0.3333333333333333 Z[(0,_0)]
   + 0.3333333333333333 Z[(1,_0)] + 0.3333333333333333 Z[(1,_1)]
   + 0.3333333333333333 Z[(2,_1)] + 0.3333333333333333 Z[(1,_2)]
   + 0.3333333333333333 Z[(2,_0)] + 0.3333333333333333 Z[(0,_2)]
   + 0.3333333333333333 Z[(2,_2)]
Subject To
 flow-constr-t0-v(0,_1): X[(0,_0),(0,_1)] + X[(1,_1),(0,_1)]
   + X[(0,_2),(0,_1)] - Y[(0,_1)] = 0
 flow-constr-t0-v(0,_0): X[(0,_1),(0,_0)] + X[(1,_0),(0,_0)] - Y[(0,_0)]
   = 0
 flow-constr-t0-v(1,_0): X[(0,_0),(1,_0)] + X[(1,_1),(1,_0)]
   + X[(2,_0),(1,_0)] - Y[(1,_0)] = 0
 flow-constr-t0-v(1,_1): X[(0,_1),(1,_1)] + X[(1,_2),(1,_1)]
   + X[(2,_1),(1,_1)] + X[(1,_0),(1,_1)] - Y[(1,_1)] = 0
 flow-constr-t0-v(2,_1): X[(2,_0),(2,_1)] + X[(1,_1),(2,_1)]
   + X[(2,_2),(2,_1)] - Y[(2,_1)] = 0
 flow-constr-t0-v(1,_2): X[(1,_1),(1,_2)] + X[(0,_2),(1,_2)]
   + X[(2,_2),(1,_2)] - Y[(1,_2)] = 0
 flow-constr-t0-v(2,_0): X[(2,_1),(2,_0)] + X[(1,_0),(2,_0)] - Y[(2,_0)]
   = 0
 flow-constr-t0-v(0,_2): X[(0,_1),(0,_2)] + X[(1,_2),(0,_2)] - Y[(0,_2)]
   = 0
 flow-constr-t0-v(2,_2): X[(1,_2),(2,_2)] + X[(2,_1),(2,_2)] - Y[(2,_2)]
   = 0
 flow-constr-t1-v(0,_1): X[(0,_0),(0,_1)] + X[(1,_1),(0,_1)]
   + X[(0,_2),(0,_1)] - Y[(0,_1)] = 0
 flow-constr-t1-v(0,_0): X[(0,_1),(0,_0)] + X[(1,_0),(0,_0)] - Y[(0,_0)]
   = 0
 flow-constr-t1-v(1,_0): X[(0,_0),(1,_0)] + X[(1,_1),(1,_0)]
   + X[(2,_0),(1,_0)] - Y[(1,_0)] = 0
 flow-constr-t1-v(1,_1): X[(0,_1),(1,_1)] + X[(1,_2),(1,_1)]
   + X[(2,_1),(1,_1)] + X[(1,_0),(1,_1)] - Y[(1,_1)] = 0
 flow-constr-t1-v(2,_1): X[(2,_0),(2,_1)] + X[(1,_1),(2,_1)]
   + X[(2,_2),(2,_1)] - Y[(2,_1)] = 0
 flow-constr-t1-v(1,_2): X[(1,_1),(1,_2)] + X[(0,_2),(1,_2)]
   + X[(2,_2),(1,_2)] - Y[(1,_2)] = 0
 flow-constr-t1-v(2,_0): X[(2,_1),(2,_0)] + X[(1,_0),(2,_0)] - Y[(2,_0)]
   = 0
 flow-constr-t1-v(0,_2): X[(0,_1),(0,_2)] + X[(1,_2),(0,_2)] - Y[(0,_2)]
   = 0
 flow-constr-t1-v(2,_2): X[(1,_2),(2,_2)] + X[(2,_1),(2,_2)] - Y[(2,_2)]
   = 0
 flow-constr-t2-v(0,_1): X[(0,_0),(0,_1)] + X[(1,_1),(0,_1)]
   + X[(0,_2),(0,_1)] - Y[(0,_1)] = 0
 flow-constr-t2-v(0,_0): X[(0,_1),(0,_0)] + X[(1,_0),(0,_0)] - Y[(0,_0)]
   = 0
 flow-constr-t2-v(1,_0): X[(0,_0),(1,_0)] + X[(1,_1),(1,_0)]
   + X[(2,_0),(1,_0)] - Y[(1,_0)] = 0
 flow-constr-t2-v(1,_1): X[(0,_1),(1,_1)] + X[(1,_2),(1,_1)]
   + X[(2,_1),(1,_1)] + X[(1,_0),(1,_1)] - Y[(1,_1)] = 0
 flow-constr-t2-v(2,_1): X[(2,_0),(2,_1)] + X[(1,_1),(2,_1)]
   + X[(2,_2),(2,_1)] - Y[(2,_1)] = 0
 flow-constr-t2-v(1,_2): X[(1,_1),(1,_2)] + X[(0,_2),(1,_2)]
   + X[(2,_2),(1,_2)] - Y[(1,_2)] = 0
 flow-constr-t2-v(2,_0): X[(2,_1),(2,_0)] + X[(1,_0),(2,_0)] - Y[(2,_0)]
   = 0
 flow-constr-t2-v(0,_2): X[(0,_1),(0,_2)] + X[(1,_2),(0,_2)] - Y[(0,_2)]
   = 0
 flow-constr-t2-v(2,_2): X[(1,_2),(2,_2)] + X[(2,_1),(2,_2)] - Y[(2,_2)]
   = 0
 root-constr-t0: Y[(0,_1)] + Y[(1,_0)] >= 1
 root-constr-t1: Y[(0,_1)] + Y[(1,_0)] >= 1
 root-constr-t2: Y[(0,_1)] + Y[(1,_0)] >= 1
 weight-constr-t2: Z[(0,_1)] + Z[(0,_1)] + Z[(0,_1)] >= 3
 weight-constr-t2: Z[(1,_0)] + Z[(1,_0)] + Z[(1,_0)] >= 3
 weight-constr-t2: Z[(1,_1)] + Z[(1,_1)] + Z[(1,_1)] >= 3
 weight-constr-t2: Z[(2,_1)] + Z[(2,_1)] + Z[(2,_1)] >= 3
 weight-constr-t2: Z[(1,_2)] + Z[(1,_2)] + Z[(1,_2)] >= 3
 weight-constr-t2: Z[(2,_0)] + Z[(2,_0)] + Z[(2,_0)] >= 3
 weight-constr-t2: Z[(0,_2)] + Z[(0,_2)] + Z[(0,_2)] >= 3
 weight-constr-t2: Z[(2,_2)] + Z[(2,_2)] + Z[(2,_2)] >= 3
 strat-constr-t0-a(1,_1)(0,_1): - 32 X[(1,_1),(0,_1)] - 2 Z[(0,_1)]
   + Z[(1,_1)] >= -32
 strat-constr-t0-a(0,_1)(1,_1): - 32 X[(0,_1),(1,_1)] + Z[(0,_1)]
   - 2 Z[(1,_1)] >= -32
 strat-constr-t0-a(0,_2)(0,_1): - 32 X[(0,_2),(0,_1)] - 2 Z[(0,_1)]
   + Z[(0,_2)] >= -32
 strat-constr-t0-a(0,_1)(0,_2): - 32 X[(0,_1),(0,_2)] + Z[(0,_1)]
   - 2 Z[(0,_2)] >= -32
 strat-constr-t0-a(1,_1)(1,_2): - 32 X[(1,_1),(1,_2)] + Z[(1,_1)]
   - 2 Z[(1,_2)] >= -32
 strat-constr-t0-a(1,_2)(1,_1): - 32 X[(1,_2),(1,_1)] - 2 Z[(1,_1)]
   + Z[(1,_2)] >= -32
 strat-constr-t0-a(0,_2)(1,_2): - 32 X[(0,_2),(1,_2)] - 2 Z[(1,_2)]
   + Z[(0,_2)] >= -32
 strat-constr-t0-a(1,_2)(0,_2): - 32 X[(1,_2),(0,_2)] + Z[(1,_2)]
   - 2 Z[(0,_2)] >= -32
 strat-constr-t0-a(2,_2)(1,_2): - 32 X[(2,_2),(1,_2)] - 2 Z[(1,_2)]
   + Z[(2,_2)] >= -32
 strat-constr-t0-a(1,_2)(2,_2): - 32 X[(1,_2),(2,_2)] + Z[(1,_2)]
   - 2 Z[(2,_2)] >= -32
 strat-constr-t0-a(2,_0)(2,_1): - 32 X[(2,_0),(2,_1)] - 2 Z[(2,_1)]
   + Z[(2,_0)] >= -32
 strat-constr-t0-a(2,_1)(2,_0): - 32 X[(2,_1),(2,_0)] + Z[(2,_1)]
   - 2 Z[(2,_0)] >= -32
 strat-constr-t0-a(1,_1)(2,_1): - 32 X[(1,_1),(2,_1)] + Z[(1,_1)]
   - 2 Z[(2,_1)] >= -32
 strat-constr-t0-a(2,_1)(1,_1): - 32 X[(2,_1),(1,_1)] - 2 Z[(1,_1)]
   + Z[(2,_1)] >= -32
 strat-constr-t0-a(2,_2)(2,_1): - 32 X[(2,_2),(2,_1)] - 2 Z[(2,_1)]
   + Z[(2,_2)] >= -32
 strat-constr-t0-a(2,_1)(2,_2): - 32 X[(2,_1),(2,_2)] + Z[(2,_1)]
   - 2 Z[(2,_2)] >= -32
 strat-constr-t0-a(1,_0)(1,_1): - 32 X[(1,_0),(1,_1)] + Z[(1,_0)]
   - 2 Z[(1,_1)] >= -32
 strat-constr-t0-a(1,_1)(1,_0): - 32 X[(1,_1),(1,_0)] - 2 Z[(1,_0)]
   + Z[(1,_1)] >= -32
 strat-constr-t0-a(1,_0)(2,_0): - 32 X[(1,_0),(2,_0)] + Z[(1,_0)]
   - 2 Z[(2,_0)] >= -32
 strat-constr-t0-a(2,_0)(1,_0): - 32 X[(2,_0),(1,_0)] - 2 Z[(1,_0)]
   + Z[(2,_0)] >= -32
 strat-constr-t1-a(1,_1)(0,_1): - 32 X[(1,_1),(0,_1)] - 2 Z[(0,_1)]
   + Z[(1,_1)] >= -32
 strat-constr-t1-a(0,_1)(1,_1): - 32 X[(0,_1),(1,_1)] + Z[(0,_1)]
   - 2 Z[(1,_1)] >= -32
 strat-constr-t1-a(0,_2)(0,_1): - 32 X[(0,_2),(0,_1)] - 2 Z[(0,_1)]
   + Z[(0,_2)] >= -32
 strat-constr-t1-a(0,_1)(0,_2): - 32 X[(0,_1),(0,_2)] + Z[(0,_1)]
   - 2 Z[(0,_2)] >= -32
 strat-constr-t1-a(1,_1)(1,_2): - 32 X[(1,_1),(1,_2)] + Z[(1,_1)]
   - 2 Z[(1,_2)] >= -32
 strat-constr-t1-a(1,_2)(1,_1): - 32 X[(1,_2),(1,_1)] - 2 Z[(1,_1)]
   + Z[(1,_2)] >= -32
 strat-constr-t1-a(0,_2)(1,_2): - 32 X[(0,_2),(1,_2)] - 2 Z[(1,_2)]
   + Z[(0,_2)] >= -32
 strat-constr-t1-a(1,_2)(0,_2): - 32 X[(1,_2),(0,_2)] + Z[(1,_2)]
   - 2 Z[(0,_2)] >= -32
 strat-constr-t1-a(2,_2)(1,_2): - 32 X[(2,_2),(1,_2)] - 2 Z[(1,_2)]
   + Z[(2,_2)] >= -32
 strat-constr-t1-a(1,_2)(2,_2): - 32 X[(1,_2),(2,_2)] + Z[(1,_2)]
   - 2 Z[(2,_2)] >= -32
 strat-constr-t1-a(2,_0)(2,_1): - 32 X[(2,_0),(2,_1)] - 2 Z[(2,_1)]
   + Z[(2,_0)] >= -32
 strat-constr-t1-a(2,_1)(2,_0): - 32 X[(2,_1),(2,_0)] + Z[(2,_1)]
   - 2 Z[(2,_0)] >= -32
 strat-constr-t1-a(1,_1)(2,_1): - 32 X[(1,_1),(2,_1)] + Z[(1,_1)]
   - 2 Z[(2,_1)] >= -32
 strat-constr-t1-a(2,_1)(1,_1): - 32 X[(2,_1),(1,_1)] - 2 Z[(1,_1)]
   + Z[(2,_1)] >= -32
 strat-constr-t1-a(2,_2)(2,_1): - 32 X[(2,_2),(2,_1)] - 2 Z[(2,_1)]
   + Z[(2,_2)] >= -32
 strat-constr-t1-a(2,_1)(2,_2): - 32 X[(2,_1),(2,_2)] + Z[(2,_1)]
   - 2 Z[(2,_2)] >= -32
 strat-constr-t1-a(1,_0)(1,_1): - 32 X[(1,_0),(1,_1)] + Z[(1,_0)]
   - 2 Z[(1,_1)] >= -32
 strat-constr-t1-a(1,_1)(1,_0): - 32 X[(1,_1),(1,_0)] - 2 Z[(1,_0)]
   + Z[(1,_1)] >= -32
 strat-constr-t1-a(1,_0)(2,_0): - 32 X[(1,_0),(2,_0)] + Z[(1,_0)]
   - 2 Z[(2,_0)] >= -32
 strat-constr-t1-a(2,_0)(1,_0): - 32 X[(2,_0),(1,_0)] - 2 Z[(1,_0)]
   + Z[(2,_0)] >= -32
 strat-constr-t2-a(1,_1)(0,_1): - 32 X[(1,_1),(0,_1)] - 2 Z[(0,_1)]
   + Z[(1,_1)] >= -32
 strat-constr-t2-a(0,_1)(1,_1): - 32 X[(0,_1),(1,_1)] + Z[(0,_1)]
   - 2 Z[(1,_1)] >= -32
 strat-constr-t2-a(0,_2)(0,_1): - 32 X[(0,_2),(0,_1)] - 2 Z[(0,_1)]
   + Z[(0,_2)] >= -32
 strat-constr-t2-a(0,_1)(0,_2): - 32 X[(0,_1),(0,_2)] + Z[(0,_1)]
   - 2 Z[(0,_2)] >= -32
 strat-constr-t2-a(1,_1)(1,_2): - 32 X[(1,_1),(1,_2)] + Z[(1,_1)]
   - 2 Z[(1,_2)] >= -32
 strat-constr-t2-a(1,_2)(1,_1): - 32 X[(1,_2),(1,_1)] - 2 Z[(1,_1)]
   + Z[(1,_2)] >= -32
 strat-constr-t2-a(0,_2)(1,_2): - 32 X[(0,_2),(1,_2)] - 2 Z[(1,_2)]
   + Z[(0,_2)] >= -32
 strat-constr-t2-a(1,_2)(0,_2): - 32 X[(1,_2),(0,_2)] + Z[(1,_2)]
   - 2 Z[(0,_2)] >= -32
 strat-constr-t2-a(2,_2)(1,_2): - 32 X[(2,_2),(1,_2)] - 2 Z[(1,_2)]
   + Z[(2,_2)] >= -32
 strat-constr-t2-a(1,_2)(2,_2): - 32 X[(1,_2),(2,_2)] + Z[(1,_2)]
   - 2 Z[(2,_2)] >= -32
 strat-constr-t2-a(2,_0)(2,_1): - 32 X[(2,_0),(2,_1)] - 2 Z[(2,_1)]
   + Z[(2,_0)] >= -32
 strat-constr-t2-a(2,_1)(2,_0): - 32 X[(2,_1),(2,_0)] + Z[(2,_1)]
   - 2 Z[(2,_0)] >= -32
 strat-constr-t2-a(1,_1)(2,_1): - 32 X[(1,_1),(2,_1)] + Z[(1,_1)]
   - 2 Z[(2,_1)] >= -32
 strat-constr-t2-a(2,_1)(1,_1): - 32 X[(2,_1),(1,_1)] - 2 Z[(1,_1)]
   + Z[(2,_1)] >= -32
 strat-constr-t2-a(2,_2)(2,_1): - 32 X[(2,_2),(2,_1)] - 2 Z[(2,_1)]
   + Z[(2,_2)] >= -32
 strat-constr-t2-a(2,_1)(2,_2): - 32 X[(2,_1),(2,_2)] + Z[(2,_1)]
   - 2 Z[(2,_2)] >= -32
 strat-constr-t2-a(1,_0)(1,_1): - 32 X[(1,_0),(1,_1)] + Z[(1,_0)]
   - 2 Z[(1,_1)] >= -32
 strat-constr-t2-a(1,_1)(1,_0): - 32 X[(1,_1),(1,_0)] - 2 Z[(1,_0)]
   + Z[(1,_1)] >= -32
 strat-constr-t2-a(1,_0)(2,_0): - 32 X[(1,_0),(2,_0)] + Z[(1,_0)]
   - 2 Z[(2,_0)] >= -32
 strat-constr-t2-a(2,_0)(1,_0): - 32 X[(2,_0),(1,_0)] - 2 Z[(1,_0)]
   + Z[(2,_0)] >= -32
 weight-constr2-t0-v(0,_1): - 16 Y[(0,_1)] + Z[(0,_1)] <= 0
 weight-constr3-t0-v(0,_1): Y[(0,_1)] - Z[(0,_1)] <= 0
 weight-constr2-t0-v(0,_0): - 16 Y[(0,_0)] + Z[(0,_0)] <= 0
 weight-constr3-t0-v(0,_0): Y[(0,_0)] - Z[(0,_0)] <= 0
 weight-constr2-t0-v(1,_0): - 16 Y[(1,_0)] + Z[(1,_0)] <= 0
 weight-constr3-t0-v(1,_0): Y[(1,_0)] - Z[(1,_0)] <= 0
 weight-constr2-t0-v(1,_1): - 16 Y[(1,_1)] + Z[(1,_1)] <= 0
 weight-constr3-t0-v(1,_1): Y[(1,_1)] - Z[(1,_1)] <= 0
 weight-constr2-t0-v(2,_1): - 16 Y[(2,_1)] + Z[(2,_1)] <= 0
 weight-constr3-t0-v(2,_1): Y[(2,_1)] - Z[(2,_1)] <= 0
 weight-constr2-t0-v(1,_2): - 16 Y[(1,_2)] + Z[(1,_2)] <= 0
 weight-constr3-t0-v(1,_2): Y[(1,_2)] - Z[(1,_2)] <= 0
 weight-constr2-t0-v(2,_0): - 16 Y[(2,_0)] + Z[(2,_0)] <= 0
 weight-constr3-t0-v(2,_0): Y[(2,_0)] - Z[(2,_0)] <= 0
 weight-constr2-t0-v(0,_2): - 16 Y[(0,_2)] + Z[(0,_2)] <= 0
 weight-constr3-t0-v(0,_2): Y[(0,_2)] - Z[(0,_2)] <= 0
 weight-constr2-t0-v(2,_2): - 16 Y[(2,_2)] + Z[(2,_2)] <= 0
 weight-constr3-t0-v(2,_2): Y[(2,_2)] - Z[(2,_2)] <= 0
 weight-constr2-t1-v(0,_1): - 16 Y[(0,_1)] + Z[(0,_1)] <= 0
 weight-constr3-t1-v(0,_1): Y[(0,_1)] - Z[(0,_1)] <= 0
 weight-constr2-t1-v(0,_0): - 16 Y[(0,_0)] + Z[(0,_0)] <= 0
 weight-constr3-t1-v(0,_0): Y[(0,_0)] - Z[(0,_0)] <= 0
 weight-constr2-t1-v(1,_0): - 16 Y[(1,_0)] + Z[(1,_0)] <= 0
 weight-constr3-t1-v(1,_0): Y[(1,_0)] - Z[(1,_0)] <= 0
 weight-constr2-t1-v(1,_1): - 16 Y[(1,_1)] + Z[(1,_1)] <= 0
 weight-constr3-t1-v(1,_1): Y[(1,_1)] - Z[(1,_1)] <= 0
 weight-constr2-t1-v(2,_1): - 16 Y[(2,_1)] + Z[(2,_1)] <= 0
 weight-constr3-t1-v(2,_1): Y[(2,_1)] - Z[(2,_1)] <= 0
 weight-constr2-t1-v(1,_2): - 16 Y[(1,_2)] + Z[(1,_2)] <= 0
 weight-constr3-t1-v(1,_2): Y[(1,_2)] - Z[(1,_2)] <= 0
 weight-constr2-t1-v(2,_0): - 16 Y[(2,_0)] + Z[(2,_0)] <= 0
 weight-constr3-t1-v(2,_0): Y[(2,_0)] - Z[(2,_0)] <= 0
 weight-constr2-t1-v(0,_2): - 16 Y[(0,_2)] + Z[(0,_2)] <= 0
 weight-constr3-t1-v(0,_2): Y[(0,_2)] - Z[(0,_2)] <= 0
 weight-constr2-t1-v(2,_2): - 16 Y[(2,_2)] + Z[(2,_2)] <= 0
 weight-constr3-t1-v(2,_2): Y[(2,_2)] - Z[(2,_2)] <= 0
 weight-constr2-t2-v(0,_1): - 16 Y[(0,_1)] + Z[(0,_1)] <= 0
 weight-constr3-t2-v(0,_1): Y[(0,_1)] - Z[(0,_1)] <= 0
 weight-constr2-t2-v(0,_0): - 16 Y[(0,_0)] + Z[(0,_0)] <= 0
 weight-constr3-t2-v(0,_0): Y[(0,_0)] - Z[(0,_0)] <= 0
 weight-constr2-t2-v(1,_0): - 16 Y[(1,_0)] + Z[(1,_0)] <= 0
 weight-constr3-t2-v(1,_0): Y[(1,_0)] - Z[(1,_0)] <= 0
 weight-constr2-t2-v(1,_1): - 16 Y[(1,_1)] + Z[(1,_1)] <= 0
 weight-constr3-t2-v(1,_1): Y[(1,_1)] - Z[(1,_1)] <= 0
 weight-constr2-t2-v(2,_1): - 16 Y[(2,_1)] + Z[(2,_1)] <= 0
 weight-constr3-t2-v(2,_1): Y[(2,_1)] - Z[(2,_1)] <= 0
 weight-constr2-t2-v(1,_2): - 16 Y[(1,_2)] + Z[(1,_2)] <= 0
 weight-constr3-t2-v(1,_2): Y[(1,_2)] - Z[(1,_2)] <= 0
 weight-constr2-t2-v(2,_0): - 16 Y[(2,_0)] + Z[(2,_0)] <= 0
 weight-constr3-t2-v(2,_0): Y[(2,_0)] - Z[(2,_0)] <= 0
 weight-constr2-t2-v(0,_2): - 16 Y[(0,_2)] + Z[(0,_2)] <= 0
 weight-constr3-t2-v(0,_2): Y[(0,_2)] - Z[(0,_2)] <= 0
 weight-constr2-t2-v(2,_2): - 16 Y[(2,_2)] + Z[(2,_2)] <= 0
 weight-constr3-t2-v(2,_2): Y[(2,_2)] - Z[(2,_2)] <= 0
Bounds
 Z[(0,_1)] <= 16
 Z[(0,_0)] <= 16
 Z[(1,_0)] <= 16
 Z[(1,_1)] <= 16
 Z[(2,_1)] <= 16
 Z[(1,_2)] <= 16
 Z[(2,_0)] <= 16
 Z[(0,_2)] <= 16
 Z[(2,_2)] <= 16
 Z[(0,_1)] <= 16
 Z[(0,_0)] <= 16
 Z[(1,_0)] <= 16
 Z[(1,_1)] <= 16
 Z[(2,_1)] <= 16
 Z[(1,_2)] <= 16
 Z[(2,_0)] <= 16
 Z[(0,_2)] <= 16
 Z[(2,_2)] <= 16
 Z[(0,_1)] <= 16
 Z[(0,_0)] <= 16
 Z[(1,_0)] <= 16
 Z[(1,_1)] <= 16
 Z[(2,_1)] <= 16
 Z[(1,_2)] <= 16
 Z[(2,_0)] <= 16
 Z[(0,_2)] <= 16
 Z[(2,_2)] <= 16
Binaries
 X[(0,_0),(0,_1)] X[(0,_1),(0,_0)] X[(1,_1),(0,_1)] X[(0,_1),(1,_1)]
 X[(0,_2),(0,_1)] X[(0,_1),(0,_2)] X[(1,_1),(1,_2)] X[(1,_2),(1,_1)]
 X[(0,_2),(1,_2)] X[(1,_2),(0,_2)] X[(2,_2),(1,_2)] X[(1,_2),(2,_2)]
 X[(1,_0),(0,_0)] X[(0,_0),(1,_0)] X[(2,_0),(2,_1)] X[(2,_1),(2,_0)]
 X[(1,_1),(2,_1)] X[(2,_1),(1,_1)] X[(2,_2),(2,_1)] X[(2,_1),(2,_2)]
 X[(1,_0),(1,_1)] X[(1,_1),(1,_0)] X[(1,_0),(2,_0)] X[(2,_0),(1,_0)]
 Y[(0,_1)] Y[(0,_0)] Y[(1,_0)] Y[(1,_1)] Y[(2,_1)] Y[(1,_2)] Y[(2,_0)]
 Y[(0,_2)] Y[(2,_2)] X[(0,_0),(0,_1)] X[(0,_1),(0,_0)] X[(1,_1),(0,_1)]
 X[(0,_1),(1,_1)] X[(0,_2),(0,_1)] X[(0,_1),(0,_2)] X[(1,_1),(1,_2)]
 X[(1,_2),(1,_1)] X[(0,_2),(1,_2)] X[(1,_2),(0,_2)] X[(2,_2),(1,_2)]
 X[(1,_2),(2,_2)] X[(1,_0),(0,_0)] X[(0,_0),(1,_0)] X[(2,_0),(2,_1)]
 X[(2,_1),(2,_0)] X[(1,_1),(2,_1)] X[(2,_1),(1,_1)] X[(2,_2),(2,_1)]
 X[(2,_1),(2,_2)] X[(1,_0),(1,_1)] X[(1,_1),(1,_0)] X[(1,_0),(2,_0)]
 X[(2,_0),(1,_0)] Y[(0,_1)] Y[(0,_0)] Y[(1,_0)] Y[(1,_1)] Y[(2,_1)]
 Y[(1,_2)] Y[(2,_0)] Y[(0,_2)] Y[(2,_2)] X[(0,_0),(0,_1)] X[(0,_1),(0,_0)]
 X[(1,_1),(0,_1)] X[(0,_1),(1,_1)] X[(0,_2),(0,_1)] X[(0,_1),(0,_2)]
 X[(1,_1),(1,_2)] X[(1,_2),(1,_1)] X[(0,_2),(1,_2)] X[(1,_2),(0,_2)]
 X[(2,_2),(1,_2)] X[(1,_2),(2,_2)] X[(1,_0),(0,_0)] X[(0,_0),(1,_0)]
 X[(2,_0),(2,_1)] X[(2,_1),(2,_0)] X[(1,_1),(2,_1)] X[(2,_1),(1,_1)]
 X[(2,_2),(2,_1)] X[(2,_1),(2,_2)] X[(1,_0),(1,_1)] X[(1,_1),(1,_0)]
 X[(1,_0),(2,_0)] X[(2,_0),(1,_0)] Y[(0,_1)] Y[(0,_0)] Y[(1,_0)] Y[(1,_1)]
 Y[(2,_1)] Y[(1,_2)] Y[(2,_0)] Y[(0,_2)] Y[(2,_2)]
Generals
 Z[(0,_1)] Z[(0,_0)] Z[(1,_0)] Z[(1,_1)] Z[(2,_1)] Z[(1,_2)] Z[(2,_0)]
 Z[(0,_2)] Z[(2,_2)] Z[(0,_1)] Z[(0,_0)] Z[(1,_0)] Z[(1,_1)] Z[(2,_1)]
 Z[(1,_2)] Z[(2,_0)] Z[(0,_2)] Z[(2,_2)] Z[(0,_1)] Z[(0,_0)] Z[(1,_0)]
 Z[(1,_1)] Z[(2,_1)] Z[(1,_2)] Z[(2,_0)] Z[(0,_2)] Z[(2,_2)]
End
