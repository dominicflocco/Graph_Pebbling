\ Model tree-optimizer
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  0.1666666666666667 Z[N11] + 0.1666666666666667 Z[N7]
   + 0.1666666666666667 Z[N8] + 0.1666666666666667 Z[N4]
   + 0.1666666666666667 Z[N2] + 0.1666666666666667 Z[N15]
   + 0.1666666666666667 Z[N5] + 0.1666666666666667 Z[N6]
   + 0.1666666666666667 Z[N3] + 0.1666666666666667 Z[N12]
   + 0.1666666666666667 Z[N14] + 0.1666666666666667 Z[N9]
   + 0.1666666666666667 Z[N10] + 0.1666666666666667 Z[N13]
   + 0.1666666666666667 Z[N1] + 0.1666666666666667 Z[N11]
   + 0.1666666666666667 Z[N7] + 0.1666666666666667 Z[N8]
   + 0.1666666666666667 Z[N4] + 0.1666666666666667 Z[N2]
   + 0.1666666666666667 Z[N15] + 0.1666666666666667 Z[N5]
   + 0.1666666666666667 Z[N6] + 0.1666666666666667 Z[N3]
   + 0.1666666666666667 Z[N12] + 0.1666666666666667 Z[N14]
   + 0.1666666666666667 Z[N9] + 0.1666666666666667 Z[N10]
   + 0.1666666666666667 Z[N13] + 0.1666666666666667 Z[N1]
   + 0.1666666666666667 Z[N11] + 0.1666666666666667 Z[N7]
   + 0.1666666666666667 Z[N8] + 0.1666666666666667 Z[N4]
   + 0.1666666666666667 Z[N2] + 0.1666666666666667 Z[N15]
   + 0.1666666666666667 Z[N5] + 0.1666666666666667 Z[N6]
   + 0.1666666666666667 Z[N3] + 0.1666666666666667 Z[N12]
   + 0.1666666666666667 Z[N14] + 0.1666666666666667 Z[N9]
   + 0.1666666666666667 Z[N10] + 0.1666666666666667 Z[N13]
   + 0.1666666666666667 Z[N1] + 0.1666666666666667 Z[N11]
   + 0.1666666666666667 Z[N7] + 0.1666666666666667 Z[N8]
   + 0.1666666666666667 Z[N4] + 0.1666666666666667 Z[N2]
   + 0.1666666666666667 Z[N15] + 0.1666666666666667 Z[N5]
   + 0.1666666666666667 Z[N6] + 0.1666666666666667 Z[N3]
   + 0.1666666666666667 Z[N12] + 0.1666666666666667 Z[N14]
   + 0.1666666666666667 Z[N9] + 0.1666666666666667 Z[N10]
   + 0.1666666666666667 Z[N13] + 0.1666666666666667 Z[N1]
   + 0.1666666666666667 Z[N11] + 0.1666666666666667 Z[N7]
   + 0.1666666666666667 Z[N8] + 0.1666666666666667 Z[N4]
   + 0.1666666666666667 Z[N2] + 0.1666666666666667 Z[N15]
   + 0.1666666666666667 Z[N5] + 0.1666666666666667 Z[N6]
   + 0.1666666666666667 Z[N3] + 0.1666666666666667 Z[N12]
   + 0.1666666666666667 Z[N14] + 0.1666666666666667 Z[N9]
   + 0.1666666666666667 Z[N10] + 0.1666666666666667 Z[N13]
   + 0.1666666666666667 Z[N1] + 0.1666666666666667 Z[N11]
   + 0.1666666666666667 Z[N7] + 0.1666666666666667 Z[N8]
   + 0.1666666666666667 Z[N4] + 0.1666666666666667 Z[N2]
   + 0.1666666666666667 Z[N15] + 0.1666666666666667 Z[N5]
   + 0.1666666666666667 Z[N6] + 0.1666666666666667 Z[N3]
   + 0.1666666666666667 Z[N12] + 0.1666666666666667 Z[N14]
   + 0.1666666666666667 Z[N9] + 0.1666666666666667 Z[N10]
   + 0.1666666666666667 Z[N13] + 0.1666666666666667 Z[N1]
Subject To
 flow-constr-t0-vN11: X[N13,N11] + X[N14,N11] - Y[N11] = 0
 flow-constr-t0-vN7: X[N1,N7] + X[N10,N7] + X[N14,N7] - Y[N7] = 0
 flow-constr-t0-vN8: X[N1,N8] + X[N3,N8] + X[N9,N8] + X[N15,N8] - Y[N8]
   = 0
 flow-constr-t0-vN4: X[N5,N4] + X[N6,N4] - Y[N4] = 0
 flow-constr-t0-vN2: X[N3,N2] + X[N5,N2] + X[N6,N2] - Y[N2] = 0
 flow-constr-t0-vN15: X[N8,N15] + X[N9,N15] - Y[N15] = 0
 flow-constr-t0-vN5: X[N2,N5] + X[N3,N5] + X[N4,N5] - Y[N5] = 0
 flow-constr-t0-vN6: X[N1,N6] + X[N2,N6] + X[N4,N6] + X[N13,N6] - Y[N6]
   = 0
 flow-constr-t0-vN3: X[N2,N3] + X[N5,N3] + X[N8,N3] - Y[N3] = 0
 flow-constr-t0-vN12: X[N9,N12] + X[N10,N12] - Y[N12] = 0
 flow-constr-t0-vN14: X[N7,N14] + X[N11,N14] + X[N13,N14] - Y[N14] = 0
 flow-constr-t0-vN9: X[N8,N9] + X[N10,N9] + X[N12,N9] + X[N15,N9] - Y[N9]
   = 0
 flow-constr-t0-vN10: X[N7,N10] + X[N9,N10] + X[N12,N10] - Y[N10] = 0
 flow-constr-t0-vN13: X[N6,N13] + X[N11,N13] + X[N14,N13] - Y[N13] = 0
 flow-constr-t0-vN1: X[N6,N1] + X[N7,N1] + X[N8,N1] - Y[N1] = 0
 flow-constr-t1-vN11: X[N13,N11] + X[N14,N11] - Y[N11] = 0
 flow-constr-t1-vN7: X[N1,N7] + X[N10,N7] + X[N14,N7] - Y[N7] = 0
 flow-constr-t1-vN8: X[N1,N8] + X[N3,N8] + X[N9,N8] + X[N15,N8] - Y[N8]
   = 0
 flow-constr-t1-vN4: X[N5,N4] + X[N6,N4] - Y[N4] = 0
 flow-constr-t1-vN2: X[N3,N2] + X[N5,N2] + X[N6,N2] - Y[N2] = 0
 flow-constr-t1-vN15: X[N8,N15] + X[N9,N15] - Y[N15] = 0
 flow-constr-t1-vN5: X[N2,N5] + X[N3,N5] + X[N4,N5] - Y[N5] = 0
 flow-constr-t1-vN6: X[N1,N6] + X[N2,N6] + X[N4,N6] + X[N13,N6] - Y[N6]
   = 0
 flow-constr-t1-vN3: X[N2,N3] + X[N5,N3] + X[N8,N3] - Y[N3] = 0
 flow-constr-t1-vN12: X[N9,N12] + X[N10,N12] - Y[N12] = 0
 flow-constr-t1-vN14: X[N7,N14] + X[N11,N14] + X[N13,N14] - Y[N14] = 0
 flow-constr-t1-vN9: X[N8,N9] + X[N10,N9] + X[N12,N9] + X[N15,N9] - Y[N9]
   = 0
 flow-constr-t1-vN10: X[N7,N10] + X[N9,N10] + X[N12,N10] - Y[N10] = 0
 flow-constr-t1-vN13: X[N6,N13] + X[N11,N13] + X[N14,N13] - Y[N13] = 0
 flow-constr-t1-vN1: X[N6,N1] + X[N7,N1] + X[N8,N1] - Y[N1] = 0
 flow-constr-t2-vN11: X[N13,N11] + X[N14,N11] - Y[N11] = 0
 flow-constr-t2-vN7: X[N1,N7] + X[N10,N7] + X[N14,N7] - Y[N7] = 0
 flow-constr-t2-vN8: X[N1,N8] + X[N3,N8] + X[N9,N8] + X[N15,N8] - Y[N8]
   = 0
 flow-constr-t2-vN4: X[N5,N4] + X[N6,N4] - Y[N4] = 0
 flow-constr-t2-vN2: X[N3,N2] + X[N5,N2] + X[N6,N2] - Y[N2] = 0
 flow-constr-t2-vN15: X[N8,N15] + X[N9,N15] - Y[N15] = 0
 flow-constr-t2-vN5: X[N2,N5] + X[N3,N5] + X[N4,N5] - Y[N5] = 0
 flow-constr-t2-vN6: X[N1,N6] + X[N2,N6] + X[N4,N6] + X[N13,N6] - Y[N6]
   = 0
 flow-constr-t2-vN3: X[N2,N3] + X[N5,N3] + X[N8,N3] - Y[N3] = 0
 flow-constr-t2-vN12: X[N9,N12] + X[N10,N12] - Y[N12] = 0
 flow-constr-t2-vN14: X[N7,N14] + X[N11,N14] + X[N13,N14] - Y[N14] = 0
 flow-constr-t2-vN9: X[N8,N9] + X[N10,N9] + X[N12,N9] + X[N15,N9] - Y[N9]
   = 0
 flow-constr-t2-vN10: X[N7,N10] + X[N9,N10] + X[N12,N10] - Y[N10] = 0
 flow-constr-t2-vN13: X[N6,N13] + X[N11,N13] + X[N14,N13] - Y[N13] = 0
 flow-constr-t2-vN1: X[N6,N1] + X[N7,N1] + X[N8,N1] - Y[N1] = 0
 flow-constr-t3-vN11: X[N13,N11] + X[N14,N11] - Y[N11] = 0
 flow-constr-t3-vN7: X[N1,N7] + X[N10,N7] + X[N14,N7] - Y[N7] = 0
 flow-constr-t3-vN8: X[N1,N8] + X[N3,N8] + X[N9,N8] + X[N15,N8] - Y[N8]
   = 0
 flow-constr-t3-vN4: X[N5,N4] + X[N6,N4] - Y[N4] = 0
 flow-constr-t3-vN2: X[N3,N2] + X[N5,N2] + X[N6,N2] - Y[N2] = 0
 flow-constr-t3-vN15: X[N8,N15] + X[N9,N15] - Y[N15] = 0
 flow-constr-t3-vN5: X[N2,N5] + X[N3,N5] + X[N4,N5] - Y[N5] = 0
 flow-constr-t3-vN6: X[N1,N6] + X[N2,N6] + X[N4,N6] + X[N13,N6] - Y[N6]
   = 0
 flow-constr-t3-vN3: X[N2,N3] + X[N5,N3] + X[N8,N3] - Y[N3] = 0
 flow-constr-t3-vN12: X[N9,N12] + X[N10,N12] - Y[N12] = 0
 flow-constr-t3-vN14: X[N7,N14] + X[N11,N14] + X[N13,N14] - Y[N14] = 0
 flow-constr-t3-vN9: X[N8,N9] + X[N10,N9] + X[N12,N9] + X[N15,N9] - Y[N9]
   = 0
 flow-constr-t3-vN10: X[N7,N10] + X[N9,N10] + X[N12,N10] - Y[N10] = 0
 flow-constr-t3-vN13: X[N6,N13] + X[N11,N13] + X[N14,N13] - Y[N13] = 0
 flow-constr-t3-vN1: X[N6,N1] + X[N7,N1] + X[N8,N1] - Y[N1] = 0
 flow-constr-t4-vN11: X[N13,N11] + X[N14,N11] - Y[N11] = 0
 flow-constr-t4-vN7: X[N1,N7] + X[N10,N7] + X[N14,N7] - Y[N7] = 0
 flow-constr-t4-vN8: X[N1,N8] + X[N3,N8] + X[N9,N8] + X[N15,N8] - Y[N8]
   = 0
 flow-constr-t4-vN4: X[N5,N4] + X[N6,N4] - Y[N4] = 0
 flow-constr-t4-vN2: X[N3,N2] + X[N5,N2] + X[N6,N2] - Y[N2] = 0
 flow-constr-t4-vN15: X[N8,N15] + X[N9,N15] - Y[N15] = 0
 flow-constr-t4-vN5: X[N2,N5] + X[N3,N5] + X[N4,N5] - Y[N5] = 0
 flow-constr-t4-vN6: X[N1,N6] + X[N2,N6] + X[N4,N6] + X[N13,N6] - Y[N6]
   = 0
 flow-constr-t4-vN3: X[N2,N3] + X[N5,N3] + X[N8,N3] - Y[N3] = 0
 flow-constr-t4-vN12: X[N9,N12] + X[N10,N12] - Y[N12] = 0
 flow-constr-t4-vN14: X[N7,N14] + X[N11,N14] + X[N13,N14] - Y[N14] = 0
 flow-constr-t4-vN9: X[N8,N9] + X[N10,N9] + X[N12,N9] + X[N15,N9] - Y[N9]
   = 0
 flow-constr-t4-vN10: X[N7,N10] + X[N9,N10] + X[N12,N10] - Y[N10] = 0
 flow-constr-t4-vN13: X[N6,N13] + X[N11,N13] + X[N14,N13] - Y[N13] = 0
 flow-constr-t4-vN1: X[N6,N1] + X[N7,N1] + X[N8,N1] - Y[N1] = 0
 flow-constr-t5-vN11: X[N13,N11] + X[N14,N11] - Y[N11] = 0
 flow-constr-t5-vN7: X[N1,N7] + X[N10,N7] + X[N14,N7] - Y[N7] = 0
 flow-constr-t5-vN8: X[N1,N8] + X[N3,N8] + X[N9,N8] + X[N15,N8] - Y[N8]
   = 0
 flow-constr-t5-vN4: X[N5,N4] + X[N6,N4] - Y[N4] = 0
 flow-constr-t5-vN2: X[N3,N2] + X[N5,N2] + X[N6,N2] - Y[N2] = 0
 flow-constr-t5-vN15: X[N8,N15] + X[N9,N15] - Y[N15] = 0
 flow-constr-t5-vN5: X[N2,N5] + X[N3,N5] + X[N4,N5] - Y[N5] = 0
 flow-constr-t5-vN6: X[N1,N6] + X[N2,N6] + X[N4,N6] + X[N13,N6] - Y[N6]
   = 0
 flow-constr-t5-vN3: X[N2,N3] + X[N5,N3] + X[N8,N3] - Y[N3] = 0
 flow-constr-t5-vN12: X[N9,N12] + X[N10,N12] - Y[N12] = 0
 flow-constr-t5-vN14: X[N7,N14] + X[N11,N14] + X[N13,N14] - Y[N14] = 0
 flow-constr-t5-vN9: X[N8,N9] + X[N10,N9] + X[N12,N9] + X[N15,N9] - Y[N9]
   = 0
 flow-constr-t5-vN10: X[N7,N10] + X[N9,N10] + X[N12,N10] - Y[N10] = 0
 flow-constr-t5-vN13: X[N6,N13] + X[N11,N13] + X[N14,N13] - Y[N13] = 0
 flow-constr-t5-vN1: X[N6,N1] + X[N7,N1] + X[N8,N1] - Y[N1] = 0
 root-constr-t0: Y[N14] + Y[N13] >= 1
 root-constr-t1: Y[N14] + Y[N13] >= 1
 root-constr-t2: Y[N14] + Y[N13] >= 1
 root-constr-t3: Y[N14] + Y[N13] >= 1
 root-constr-t4: Y[N14] + Y[N13] >= 1
 root-constr-t5: Y[N14] + Y[N13] >= 1
 weight-constr-t5: Z[N7] + Z[N7] + Z[N7] + Z[N7] + Z[N7] + Z[N7] >= 6
 weight-constr-t5: Z[N8] + Z[N8] + Z[N8] + Z[N8] + Z[N8] + Z[N8] >= 6
 weight-constr-t5: Z[N4] + Z[N4] + Z[N4] + Z[N4] + Z[N4] + Z[N4] >= 6
 weight-constr-t5: Z[N2] + Z[N2] + Z[N2] + Z[N2] + Z[N2] + Z[N2] >= 6
 weight-constr-t5: Z[N15] + Z[N15] + Z[N15] + Z[N15] + Z[N15] + Z[N15]
   >= 6
 weight-constr-t5: Z[N5] + Z[N5] + Z[N5] + Z[N5] + Z[N5] + Z[N5] >= 6
 weight-constr-t5: Z[N6] + Z[N6] + Z[N6] + Z[N6] + Z[N6] + Z[N6] >= 6
 weight-constr-t5: Z[N3] + Z[N3] + Z[N3] + Z[N3] + Z[N3] + Z[N3] >= 6
 weight-constr-t5: Z[N12] + Z[N12] + Z[N12] + Z[N12] + Z[N12] + Z[N12]
   >= 6
 weight-constr-t5: Z[N14] + Z[N14] + Z[N14] + Z[N14] + Z[N14] + Z[N14]
   >= 6
 weight-constr-t5: Z[N9] + Z[N9] + Z[N9] + Z[N9] + Z[N9] + Z[N9] >= 6
 weight-constr-t5: Z[N10] + Z[N10] + Z[N10] + Z[N10] + Z[N10] + Z[N10]
   >= 6
 weight-constr-t5: Z[N13] + Z[N13] + Z[N13] + Z[N13] + Z[N13] + Z[N13]
   >= 6
 weight-constr-t5: Z[N1] + Z[N1] + Z[N1] + Z[N1] + Z[N1] + Z[N1] >= 6
 strat-constr-t0-aN6N1: - 1.048576e+06 X[N6,N1] + Z[N6] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t0-aN1N6: - 1.048576e+06 X[N1,N6] - 2 Z[N6] + Z[N1]
   >= -1.048576e+06
 strat-constr-t0-aN7N1: - 1.048576e+06 X[N7,N1] + Z[N7] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t0-aN1N7: - 1.048576e+06 X[N1,N7] - 2 Z[N7] + Z[N1]
   >= -1.048576e+06
 strat-constr-t0-aN8N1: - 1.048576e+06 X[N8,N1] + Z[N8] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t0-aN1N8: - 1.048576e+06 X[N1,N8] - 2 Z[N8] + Z[N1]
   >= -1.048576e+06
 strat-constr-t0-aN3N2: - 1.048576e+06 X[N3,N2] - 2 Z[N2] + Z[N3]
   >= -1.048576e+06
 strat-constr-t0-aN2N3: - 1.048576e+06 X[N2,N3] + Z[N2] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t0-aN5N2: - 1.048576e+06 X[N5,N2] - 2 Z[N2] + Z[N5]
   >= -1.048576e+06
 strat-constr-t0-aN2N5: - 1.048576e+06 X[N2,N5] + Z[N2] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t0-aN6N2: - 1.048576e+06 X[N6,N2] - 2 Z[N2] + Z[N6]
   >= -1.048576e+06
 strat-constr-t0-aN2N6: - 1.048576e+06 X[N2,N6] + Z[N2] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t0-aN5N3: - 1.048576e+06 X[N5,N3] + Z[N5] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t0-aN3N5: - 1.048576e+06 X[N3,N5] - 2 Z[N5] + Z[N3]
   >= -1.048576e+06
 strat-constr-t0-aN8N3: - 1.048576e+06 X[N8,N3] + Z[N8] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t0-aN3N8: - 1.048576e+06 X[N3,N8] - 2 Z[N8] + Z[N3]
   >= -1.048576e+06
 strat-constr-t0-aN5N4: - 1.048576e+06 X[N5,N4] - 2 Z[N4] + Z[N5]
   >= -1.048576e+06
 strat-constr-t0-aN4N5: - 1.048576e+06 X[N4,N5] + Z[N4] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t0-aN6N4: - 1.048576e+06 X[N6,N4] - 2 Z[N4] + Z[N6]
   >= -1.048576e+06
 strat-constr-t0-aN4N6: - 1.048576e+06 X[N4,N6] + Z[N4] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t0-aN13N6: - 1.048576e+06 X[N13,N6] - 2 Z[N6] + Z[N13]
   >= -1.048576e+06
 strat-constr-t0-aN6N13: - 1.048576e+06 X[N6,N13] + Z[N6] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t0-aN10N7: - 1.048576e+06 X[N10,N7] - 2 Z[N7] + Z[N10]
   >= -1.048576e+06
 strat-constr-t0-aN7N10: - 1.048576e+06 X[N7,N10] + Z[N7] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t0-aN14N7: - 1.048576e+06 X[N14,N7] - 2 Z[N7] + Z[N14]
   >= -1.048576e+06
 strat-constr-t0-aN7N14: - 1.048576e+06 X[N7,N14] + Z[N7] - 2 Z[N14]
   >= -1.048576e+06
 strat-constr-t0-aN9N8: - 1.048576e+06 X[N9,N8] - 2 Z[N8] + Z[N9]
   >= -1.048576e+06
 strat-constr-t0-aN8N9: - 1.048576e+06 X[N8,N9] + Z[N8] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t0-aN15N8: - 1.048576e+06 X[N15,N8] - 2 Z[N8] + Z[N15]
   >= -1.048576e+06
 strat-constr-t0-aN8N15: - 1.048576e+06 X[N8,N15] + Z[N8] - 2 Z[N15]
   >= -1.048576e+06
 strat-constr-t0-aN10N9: - 1.048576e+06 X[N10,N9] - 2 Z[N9] + Z[N10]
   >= -1.048576e+06
 strat-constr-t0-aN9N10: - 1.048576e+06 X[N9,N10] + Z[N9] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t0-aN12N9: - 1.048576e+06 X[N12,N9] + Z[N12] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t0-aN9N12: - 1.048576e+06 X[N9,N12] - 2 Z[N12] + Z[N9]
   >= -1.048576e+06
 strat-constr-t0-aN15N9: - 1.048576e+06 X[N15,N9] + Z[N15] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t0-aN9N15: - 1.048576e+06 X[N9,N15] - 2 Z[N15] + Z[N9]
   >= -1.048576e+06
 strat-constr-t0-aN12N10: - 1.048576e+06 X[N12,N10] + Z[N12] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t0-aN10N12: - 1.048576e+06 X[N10,N12] - 2 Z[N12] + Z[N10]
   >= -1.048576e+06
 strat-constr-t0-aN14N13: - 1.048576e+06 X[N14,N13] + Z[N14] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t0-aN13N14: - 1.048576e+06 X[N13,N14] - 2 Z[N14] + Z[N13]
   >= -1.048576e+06
 strat-constr-t1-aN6N1: - 1.048576e+06 X[N6,N1] + Z[N6] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t1-aN1N6: - 1.048576e+06 X[N1,N6] - 2 Z[N6] + Z[N1]
   >= -1.048576e+06
 strat-constr-t1-aN7N1: - 1.048576e+06 X[N7,N1] + Z[N7] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t1-aN1N7: - 1.048576e+06 X[N1,N7] - 2 Z[N7] + Z[N1]
   >= -1.048576e+06
 strat-constr-t1-aN8N1: - 1.048576e+06 X[N8,N1] + Z[N8] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t1-aN1N8: - 1.048576e+06 X[N1,N8] - 2 Z[N8] + Z[N1]
   >= -1.048576e+06
 strat-constr-t1-aN3N2: - 1.048576e+06 X[N3,N2] - 2 Z[N2] + Z[N3]
   >= -1.048576e+06
 strat-constr-t1-aN2N3: - 1.048576e+06 X[N2,N3] + Z[N2] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t1-aN5N2: - 1.048576e+06 X[N5,N2] - 2 Z[N2] + Z[N5]
   >= -1.048576e+06
 strat-constr-t1-aN2N5: - 1.048576e+06 X[N2,N5] + Z[N2] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t1-aN6N2: - 1.048576e+06 X[N6,N2] - 2 Z[N2] + Z[N6]
   >= -1.048576e+06
 strat-constr-t1-aN2N6: - 1.048576e+06 X[N2,N6] + Z[N2] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t1-aN5N3: - 1.048576e+06 X[N5,N3] + Z[N5] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t1-aN3N5: - 1.048576e+06 X[N3,N5] - 2 Z[N5] + Z[N3]
   >= -1.048576e+06
 strat-constr-t1-aN8N3: - 1.048576e+06 X[N8,N3] + Z[N8] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t1-aN3N8: - 1.048576e+06 X[N3,N8] - 2 Z[N8] + Z[N3]
   >= -1.048576e+06
 strat-constr-t1-aN5N4: - 1.048576e+06 X[N5,N4] - 2 Z[N4] + Z[N5]
   >= -1.048576e+06
 strat-constr-t1-aN4N5: - 1.048576e+06 X[N4,N5] + Z[N4] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t1-aN6N4: - 1.048576e+06 X[N6,N4] - 2 Z[N4] + Z[N6]
   >= -1.048576e+06
 strat-constr-t1-aN4N6: - 1.048576e+06 X[N4,N6] + Z[N4] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t1-aN13N6: - 1.048576e+06 X[N13,N6] - 2 Z[N6] + Z[N13]
   >= -1.048576e+06
 strat-constr-t1-aN6N13: - 1.048576e+06 X[N6,N13] + Z[N6] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t1-aN10N7: - 1.048576e+06 X[N10,N7] - 2 Z[N7] + Z[N10]
   >= -1.048576e+06
 strat-constr-t1-aN7N10: - 1.048576e+06 X[N7,N10] + Z[N7] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t1-aN14N7: - 1.048576e+06 X[N14,N7] - 2 Z[N7] + Z[N14]
   >= -1.048576e+06
 strat-constr-t1-aN7N14: - 1.048576e+06 X[N7,N14] + Z[N7] - 2 Z[N14]
   >= -1.048576e+06
 strat-constr-t1-aN9N8: - 1.048576e+06 X[N9,N8] - 2 Z[N8] + Z[N9]
   >= -1.048576e+06
 strat-constr-t1-aN8N9: - 1.048576e+06 X[N8,N9] + Z[N8] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t1-aN15N8: - 1.048576e+06 X[N15,N8] - 2 Z[N8] + Z[N15]
   >= -1.048576e+06
 strat-constr-t1-aN8N15: - 1.048576e+06 X[N8,N15] + Z[N8] - 2 Z[N15]
   >= -1.048576e+06
 strat-constr-t1-aN10N9: - 1.048576e+06 X[N10,N9] - 2 Z[N9] + Z[N10]
   >= -1.048576e+06
 strat-constr-t1-aN9N10: - 1.048576e+06 X[N9,N10] + Z[N9] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t1-aN12N9: - 1.048576e+06 X[N12,N9] + Z[N12] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t1-aN9N12: - 1.048576e+06 X[N9,N12] - 2 Z[N12] + Z[N9]
   >= -1.048576e+06
 strat-constr-t1-aN15N9: - 1.048576e+06 X[N15,N9] + Z[N15] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t1-aN9N15: - 1.048576e+06 X[N9,N15] - 2 Z[N15] + Z[N9]
   >= -1.048576e+06
 strat-constr-t1-aN12N10: - 1.048576e+06 X[N12,N10] + Z[N12] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t1-aN10N12: - 1.048576e+06 X[N10,N12] - 2 Z[N12] + Z[N10]
   >= -1.048576e+06
 strat-constr-t1-aN14N13: - 1.048576e+06 X[N14,N13] + Z[N14] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t1-aN13N14: - 1.048576e+06 X[N13,N14] - 2 Z[N14] + Z[N13]
   >= -1.048576e+06
 strat-constr-t2-aN6N1: - 1.048576e+06 X[N6,N1] + Z[N6] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t2-aN1N6: - 1.048576e+06 X[N1,N6] - 2 Z[N6] + Z[N1]
   >= -1.048576e+06
 strat-constr-t2-aN7N1: - 1.048576e+06 X[N7,N1] + Z[N7] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t2-aN1N7: - 1.048576e+06 X[N1,N7] - 2 Z[N7] + Z[N1]
   >= -1.048576e+06
 strat-constr-t2-aN8N1: - 1.048576e+06 X[N8,N1] + Z[N8] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t2-aN1N8: - 1.048576e+06 X[N1,N8] - 2 Z[N8] + Z[N1]
   >= -1.048576e+06
 strat-constr-t2-aN3N2: - 1.048576e+06 X[N3,N2] - 2 Z[N2] + Z[N3]
   >= -1.048576e+06
 strat-constr-t2-aN2N3: - 1.048576e+06 X[N2,N3] + Z[N2] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t2-aN5N2: - 1.048576e+06 X[N5,N2] - 2 Z[N2] + Z[N5]
   >= -1.048576e+06
 strat-constr-t2-aN2N5: - 1.048576e+06 X[N2,N5] + Z[N2] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t2-aN6N2: - 1.048576e+06 X[N6,N2] - 2 Z[N2] + Z[N6]
   >= -1.048576e+06
 strat-constr-t2-aN2N6: - 1.048576e+06 X[N2,N6] + Z[N2] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t2-aN5N3: - 1.048576e+06 X[N5,N3] + Z[N5] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t2-aN3N5: - 1.048576e+06 X[N3,N5] - 2 Z[N5] + Z[N3]
   >= -1.048576e+06
 strat-constr-t2-aN8N3: - 1.048576e+06 X[N8,N3] + Z[N8] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t2-aN3N8: - 1.048576e+06 X[N3,N8] - 2 Z[N8] + Z[N3]
   >= -1.048576e+06
 strat-constr-t2-aN5N4: - 1.048576e+06 X[N5,N4] - 2 Z[N4] + Z[N5]
   >= -1.048576e+06
 strat-constr-t2-aN4N5: - 1.048576e+06 X[N4,N5] + Z[N4] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t2-aN6N4: - 1.048576e+06 X[N6,N4] - 2 Z[N4] + Z[N6]
   >= -1.048576e+06
 strat-constr-t2-aN4N6: - 1.048576e+06 X[N4,N6] + Z[N4] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t2-aN13N6: - 1.048576e+06 X[N13,N6] - 2 Z[N6] + Z[N13]
   >= -1.048576e+06
 strat-constr-t2-aN6N13: - 1.048576e+06 X[N6,N13] + Z[N6] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t2-aN10N7: - 1.048576e+06 X[N10,N7] - 2 Z[N7] + Z[N10]
   >= -1.048576e+06
 strat-constr-t2-aN7N10: - 1.048576e+06 X[N7,N10] + Z[N7] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t2-aN14N7: - 1.048576e+06 X[N14,N7] - 2 Z[N7] + Z[N14]
   >= -1.048576e+06
 strat-constr-t2-aN7N14: - 1.048576e+06 X[N7,N14] + Z[N7] - 2 Z[N14]
   >= -1.048576e+06
 strat-constr-t2-aN9N8: - 1.048576e+06 X[N9,N8] - 2 Z[N8] + Z[N9]
   >= -1.048576e+06
 strat-constr-t2-aN8N9: - 1.048576e+06 X[N8,N9] + Z[N8] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t2-aN15N8: - 1.048576e+06 X[N15,N8] - 2 Z[N8] + Z[N15]
   >= -1.048576e+06
 strat-constr-t2-aN8N15: - 1.048576e+06 X[N8,N15] + Z[N8] - 2 Z[N15]
   >= -1.048576e+06
 strat-constr-t2-aN10N9: - 1.048576e+06 X[N10,N9] - 2 Z[N9] + Z[N10]
   >= -1.048576e+06
 strat-constr-t2-aN9N10: - 1.048576e+06 X[N9,N10] + Z[N9] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t2-aN12N9: - 1.048576e+06 X[N12,N9] + Z[N12] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t2-aN9N12: - 1.048576e+06 X[N9,N12] - 2 Z[N12] + Z[N9]
   >= -1.048576e+06
 strat-constr-t2-aN15N9: - 1.048576e+06 X[N15,N9] + Z[N15] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t2-aN9N15: - 1.048576e+06 X[N9,N15] - 2 Z[N15] + Z[N9]
   >= -1.048576e+06
 strat-constr-t2-aN12N10: - 1.048576e+06 X[N12,N10] + Z[N12] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t2-aN10N12: - 1.048576e+06 X[N10,N12] - 2 Z[N12] + Z[N10]
   >= -1.048576e+06
 strat-constr-t2-aN14N13: - 1.048576e+06 X[N14,N13] + Z[N14] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t2-aN13N14: - 1.048576e+06 X[N13,N14] - 2 Z[N14] + Z[N13]
   >= -1.048576e+06
 strat-constr-t3-aN6N1: - 1.048576e+06 X[N6,N1] + Z[N6] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t3-aN1N6: - 1.048576e+06 X[N1,N6] - 2 Z[N6] + Z[N1]
   >= -1.048576e+06
 strat-constr-t3-aN7N1: - 1.048576e+06 X[N7,N1] + Z[N7] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t3-aN1N7: - 1.048576e+06 X[N1,N7] - 2 Z[N7] + Z[N1]
   >= -1.048576e+06
 strat-constr-t3-aN8N1: - 1.048576e+06 X[N8,N1] + Z[N8] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t3-aN1N8: - 1.048576e+06 X[N1,N8] - 2 Z[N8] + Z[N1]
   >= -1.048576e+06
 strat-constr-t3-aN3N2: - 1.048576e+06 X[N3,N2] - 2 Z[N2] + Z[N3]
   >= -1.048576e+06
 strat-constr-t3-aN2N3: - 1.048576e+06 X[N2,N3] + Z[N2] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t3-aN5N2: - 1.048576e+06 X[N5,N2] - 2 Z[N2] + Z[N5]
   >= -1.048576e+06
 strat-constr-t3-aN2N5: - 1.048576e+06 X[N2,N5] + Z[N2] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t3-aN6N2: - 1.048576e+06 X[N6,N2] - 2 Z[N2] + Z[N6]
   >= -1.048576e+06
 strat-constr-t3-aN2N6: - 1.048576e+06 X[N2,N6] + Z[N2] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t3-aN5N3: - 1.048576e+06 X[N5,N3] + Z[N5] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t3-aN3N5: - 1.048576e+06 X[N3,N5] - 2 Z[N5] + Z[N3]
   >= -1.048576e+06
 strat-constr-t3-aN8N3: - 1.048576e+06 X[N8,N3] + Z[N8] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t3-aN3N8: - 1.048576e+06 X[N3,N8] - 2 Z[N8] + Z[N3]
   >= -1.048576e+06
 strat-constr-t3-aN5N4: - 1.048576e+06 X[N5,N4] - 2 Z[N4] + Z[N5]
   >= -1.048576e+06
 strat-constr-t3-aN4N5: - 1.048576e+06 X[N4,N5] + Z[N4] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t3-aN6N4: - 1.048576e+06 X[N6,N4] - 2 Z[N4] + Z[N6]
   >= -1.048576e+06
 strat-constr-t3-aN4N6: - 1.048576e+06 X[N4,N6] + Z[N4] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t3-aN13N6: - 1.048576e+06 X[N13,N6] - 2 Z[N6] + Z[N13]
   >= -1.048576e+06
 strat-constr-t3-aN6N13: - 1.048576e+06 X[N6,N13] + Z[N6] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t3-aN10N7: - 1.048576e+06 X[N10,N7] - 2 Z[N7] + Z[N10]
   >= -1.048576e+06
 strat-constr-t3-aN7N10: - 1.048576e+06 X[N7,N10] + Z[N7] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t3-aN14N7: - 1.048576e+06 X[N14,N7] - 2 Z[N7] + Z[N14]
   >= -1.048576e+06
 strat-constr-t3-aN7N14: - 1.048576e+06 X[N7,N14] + Z[N7] - 2 Z[N14]
   >= -1.048576e+06
 strat-constr-t3-aN9N8: - 1.048576e+06 X[N9,N8] - 2 Z[N8] + Z[N9]
   >= -1.048576e+06
 strat-constr-t3-aN8N9: - 1.048576e+06 X[N8,N9] + Z[N8] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t3-aN15N8: - 1.048576e+06 X[N15,N8] - 2 Z[N8] + Z[N15]
   >= -1.048576e+06
 strat-constr-t3-aN8N15: - 1.048576e+06 X[N8,N15] + Z[N8] - 2 Z[N15]
   >= -1.048576e+06
 strat-constr-t3-aN10N9: - 1.048576e+06 X[N10,N9] - 2 Z[N9] + Z[N10]
   >= -1.048576e+06
 strat-constr-t3-aN9N10: - 1.048576e+06 X[N9,N10] + Z[N9] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t3-aN12N9: - 1.048576e+06 X[N12,N9] + Z[N12] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t3-aN9N12: - 1.048576e+06 X[N9,N12] - 2 Z[N12] + Z[N9]
   >= -1.048576e+06
 strat-constr-t3-aN15N9: - 1.048576e+06 X[N15,N9] + Z[N15] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t3-aN9N15: - 1.048576e+06 X[N9,N15] - 2 Z[N15] + Z[N9]
   >= -1.048576e+06
 strat-constr-t3-aN12N10: - 1.048576e+06 X[N12,N10] + Z[N12] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t3-aN10N12: - 1.048576e+06 X[N10,N12] - 2 Z[N12] + Z[N10]
   >= -1.048576e+06
 strat-constr-t3-aN14N13: - 1.048576e+06 X[N14,N13] + Z[N14] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t3-aN13N14: - 1.048576e+06 X[N13,N14] - 2 Z[N14] + Z[N13]
   >= -1.048576e+06
 strat-constr-t4-aN6N1: - 1.048576e+06 X[N6,N1] + Z[N6] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t4-aN1N6: - 1.048576e+06 X[N1,N6] - 2 Z[N6] + Z[N1]
   >= -1.048576e+06
 strat-constr-t4-aN7N1: - 1.048576e+06 X[N7,N1] + Z[N7] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t4-aN1N7: - 1.048576e+06 X[N1,N7] - 2 Z[N7] + Z[N1]
   >= -1.048576e+06
 strat-constr-t4-aN8N1: - 1.048576e+06 X[N8,N1] + Z[N8] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t4-aN1N8: - 1.048576e+06 X[N1,N8] - 2 Z[N8] + Z[N1]
   >= -1.048576e+06
 strat-constr-t4-aN3N2: - 1.048576e+06 X[N3,N2] - 2 Z[N2] + Z[N3]
   >= -1.048576e+06
 strat-constr-t4-aN2N3: - 1.048576e+06 X[N2,N3] + Z[N2] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t4-aN5N2: - 1.048576e+06 X[N5,N2] - 2 Z[N2] + Z[N5]
   >= -1.048576e+06
 strat-constr-t4-aN2N5: - 1.048576e+06 X[N2,N5] + Z[N2] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t4-aN6N2: - 1.048576e+06 X[N6,N2] - 2 Z[N2] + Z[N6]
   >= -1.048576e+06
 strat-constr-t4-aN2N6: - 1.048576e+06 X[N2,N6] + Z[N2] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t4-aN5N3: - 1.048576e+06 X[N5,N3] + Z[N5] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t4-aN3N5: - 1.048576e+06 X[N3,N5] - 2 Z[N5] + Z[N3]
   >= -1.048576e+06
 strat-constr-t4-aN8N3: - 1.048576e+06 X[N8,N3] + Z[N8] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t4-aN3N8: - 1.048576e+06 X[N3,N8] - 2 Z[N8] + Z[N3]
   >= -1.048576e+06
 strat-constr-t4-aN5N4: - 1.048576e+06 X[N5,N4] - 2 Z[N4] + Z[N5]
   >= -1.048576e+06
 strat-constr-t4-aN4N5: - 1.048576e+06 X[N4,N5] + Z[N4] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t4-aN6N4: - 1.048576e+06 X[N6,N4] - 2 Z[N4] + Z[N6]
   >= -1.048576e+06
 strat-constr-t4-aN4N6: - 1.048576e+06 X[N4,N6] + Z[N4] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t4-aN13N6: - 1.048576e+06 X[N13,N6] - 2 Z[N6] + Z[N13]
   >= -1.048576e+06
 strat-constr-t4-aN6N13: - 1.048576e+06 X[N6,N13] + Z[N6] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t4-aN10N7: - 1.048576e+06 X[N10,N7] - 2 Z[N7] + Z[N10]
   >= -1.048576e+06
 strat-constr-t4-aN7N10: - 1.048576e+06 X[N7,N10] + Z[N7] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t4-aN14N7: - 1.048576e+06 X[N14,N7] - 2 Z[N7] + Z[N14]
   >= -1.048576e+06
 strat-constr-t4-aN7N14: - 1.048576e+06 X[N7,N14] + Z[N7] - 2 Z[N14]
   >= -1.048576e+06
 strat-constr-t4-aN9N8: - 1.048576e+06 X[N9,N8] - 2 Z[N8] + Z[N9]
   >= -1.048576e+06
 strat-constr-t4-aN8N9: - 1.048576e+06 X[N8,N9] + Z[N8] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t4-aN15N8: - 1.048576e+06 X[N15,N8] - 2 Z[N8] + Z[N15]
   >= -1.048576e+06
 strat-constr-t4-aN8N15: - 1.048576e+06 X[N8,N15] + Z[N8] - 2 Z[N15]
   >= -1.048576e+06
 strat-constr-t4-aN10N9: - 1.048576e+06 X[N10,N9] - 2 Z[N9] + Z[N10]
   >= -1.048576e+06
 strat-constr-t4-aN9N10: - 1.048576e+06 X[N9,N10] + Z[N9] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t4-aN12N9: - 1.048576e+06 X[N12,N9] + Z[N12] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t4-aN9N12: - 1.048576e+06 X[N9,N12] - 2 Z[N12] + Z[N9]
   >= -1.048576e+06
 strat-constr-t4-aN15N9: - 1.048576e+06 X[N15,N9] + Z[N15] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t4-aN9N15: - 1.048576e+06 X[N9,N15] - 2 Z[N15] + Z[N9]
   >= -1.048576e+06
 strat-constr-t4-aN12N10: - 1.048576e+06 X[N12,N10] + Z[N12] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t4-aN10N12: - 1.048576e+06 X[N10,N12] - 2 Z[N12] + Z[N10]
   >= -1.048576e+06
 strat-constr-t4-aN14N13: - 1.048576e+06 X[N14,N13] + Z[N14] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t4-aN13N14: - 1.048576e+06 X[N13,N14] - 2 Z[N14] + Z[N13]
   >= -1.048576e+06
 strat-constr-t5-aN6N1: - 1.048576e+06 X[N6,N1] + Z[N6] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t5-aN1N6: - 1.048576e+06 X[N1,N6] - 2 Z[N6] + Z[N1]
   >= -1.048576e+06
 strat-constr-t5-aN7N1: - 1.048576e+06 X[N7,N1] + Z[N7] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t5-aN1N7: - 1.048576e+06 X[N1,N7] - 2 Z[N7] + Z[N1]
   >= -1.048576e+06
 strat-constr-t5-aN8N1: - 1.048576e+06 X[N8,N1] + Z[N8] - 2 Z[N1]
   >= -1.048576e+06
 strat-constr-t5-aN1N8: - 1.048576e+06 X[N1,N8] - 2 Z[N8] + Z[N1]
   >= -1.048576e+06
 strat-constr-t5-aN3N2: - 1.048576e+06 X[N3,N2] - 2 Z[N2] + Z[N3]
   >= -1.048576e+06
 strat-constr-t5-aN2N3: - 1.048576e+06 X[N2,N3] + Z[N2] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t5-aN5N2: - 1.048576e+06 X[N5,N2] - 2 Z[N2] + Z[N5]
   >= -1.048576e+06
 strat-constr-t5-aN2N5: - 1.048576e+06 X[N2,N5] + Z[N2] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t5-aN6N2: - 1.048576e+06 X[N6,N2] - 2 Z[N2] + Z[N6]
   >= -1.048576e+06
 strat-constr-t5-aN2N6: - 1.048576e+06 X[N2,N6] + Z[N2] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t5-aN5N3: - 1.048576e+06 X[N5,N3] + Z[N5] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t5-aN3N5: - 1.048576e+06 X[N3,N5] - 2 Z[N5] + Z[N3]
   >= -1.048576e+06
 strat-constr-t5-aN8N3: - 1.048576e+06 X[N8,N3] + Z[N8] - 2 Z[N3]
   >= -1.048576e+06
 strat-constr-t5-aN3N8: - 1.048576e+06 X[N3,N8] - 2 Z[N8] + Z[N3]
   >= -1.048576e+06
 strat-constr-t5-aN5N4: - 1.048576e+06 X[N5,N4] - 2 Z[N4] + Z[N5]
   >= -1.048576e+06
 strat-constr-t5-aN4N5: - 1.048576e+06 X[N4,N5] + Z[N4] - 2 Z[N5]
   >= -1.048576e+06
 strat-constr-t5-aN6N4: - 1.048576e+06 X[N6,N4] - 2 Z[N4] + Z[N6]
   >= -1.048576e+06
 strat-constr-t5-aN4N6: - 1.048576e+06 X[N4,N6] + Z[N4] - 2 Z[N6]
   >= -1.048576e+06
 strat-constr-t5-aN13N6: - 1.048576e+06 X[N13,N6] - 2 Z[N6] + Z[N13]
   >= -1.048576e+06
 strat-constr-t5-aN6N13: - 1.048576e+06 X[N6,N13] + Z[N6] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t5-aN10N7: - 1.048576e+06 X[N10,N7] - 2 Z[N7] + Z[N10]
   >= -1.048576e+06
 strat-constr-t5-aN7N10: - 1.048576e+06 X[N7,N10] + Z[N7] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t5-aN14N7: - 1.048576e+06 X[N14,N7] - 2 Z[N7] + Z[N14]
   >= -1.048576e+06
 strat-constr-t5-aN7N14: - 1.048576e+06 X[N7,N14] + Z[N7] - 2 Z[N14]
   >= -1.048576e+06
 strat-constr-t5-aN9N8: - 1.048576e+06 X[N9,N8] - 2 Z[N8] + Z[N9]
   >= -1.048576e+06
 strat-constr-t5-aN8N9: - 1.048576e+06 X[N8,N9] + Z[N8] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t5-aN15N8: - 1.048576e+06 X[N15,N8] - 2 Z[N8] + Z[N15]
   >= -1.048576e+06
 strat-constr-t5-aN8N15: - 1.048576e+06 X[N8,N15] + Z[N8] - 2 Z[N15]
   >= -1.048576e+06
 strat-constr-t5-aN10N9: - 1.048576e+06 X[N10,N9] - 2 Z[N9] + Z[N10]
   >= -1.048576e+06
 strat-constr-t5-aN9N10: - 1.048576e+06 X[N9,N10] + Z[N9] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t5-aN12N9: - 1.048576e+06 X[N12,N9] + Z[N12] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t5-aN9N12: - 1.048576e+06 X[N9,N12] - 2 Z[N12] + Z[N9]
   >= -1.048576e+06
 strat-constr-t5-aN15N9: - 1.048576e+06 X[N15,N9] + Z[N15] - 2 Z[N9]
   >= -1.048576e+06
 strat-constr-t5-aN9N15: - 1.048576e+06 X[N9,N15] - 2 Z[N15] + Z[N9]
   >= -1.048576e+06
 strat-constr-t5-aN12N10: - 1.048576e+06 X[N12,N10] + Z[N12] - 2 Z[N10]
   >= -1.048576e+06
 strat-constr-t5-aN10N12: - 1.048576e+06 X[N10,N12] - 2 Z[N12] + Z[N10]
   >= -1.048576e+06
 strat-constr-t5-aN14N13: - 1.048576e+06 X[N14,N13] + Z[N14] - 2 Z[N13]
   >= -1.048576e+06
 strat-constr-t5-aN13N14: - 1.048576e+06 X[N13,N14] - 2 Z[N14] + Z[N13]
   >= -1.048576e+06
 weight-constr2-t0-vN11: - 524288 Y[N11] + Z[N11] <= 0
 weight-constr3-t0-vN11: Y[N11] - Z[N11] <= 0
 weight-constr2-t0-vN7: - 524288 Y[N7] + Z[N7] <= 0
 weight-constr3-t0-vN7: Y[N7] - Z[N7] <= 0
 weight-constr2-t0-vN8: - 524288 Y[N8] + Z[N8] <= 0
 weight-constr3-t0-vN8: Y[N8] - Z[N8] <= 0
 weight-constr2-t0-vN4: - 524288 Y[N4] + Z[N4] <= 0
 weight-constr3-t0-vN4: Y[N4] - Z[N4] <= 0
 weight-constr2-t0-vN2: - 524288 Y[N2] + Z[N2] <= 0
 weight-constr3-t0-vN2: Y[N2] - Z[N2] <= 0
 weight-constr2-t0-vN15: - 524288 Y[N15] + Z[N15] <= 0
 weight-constr3-t0-vN15: Y[N15] - Z[N15] <= 0
 weight-constr2-t0-vN5: - 524288 Y[N5] + Z[N5] <= 0
 weight-constr3-t0-vN5: Y[N5] - Z[N5] <= 0
 weight-constr2-t0-vN6: - 524288 Y[N6] + Z[N6] <= 0
 weight-constr3-t0-vN6: Y[N6] - Z[N6] <= 0
 weight-constr2-t0-vN3: - 524288 Y[N3] + Z[N3] <= 0
 weight-constr3-t0-vN3: Y[N3] - Z[N3] <= 0
 weight-constr2-t0-vN12: - 524288 Y[N12] + Z[N12] <= 0
 weight-constr3-t0-vN12: Y[N12] - Z[N12] <= 0
 weight-constr2-t0-vN14: - 524288 Y[N14] + Z[N14] <= 0
 weight-constr3-t0-vN14: Y[N14] - Z[N14] <= 0
 weight-constr2-t0-vN9: - 524288 Y[N9] + Z[N9] <= 0
 weight-constr3-t0-vN9: Y[N9] - Z[N9] <= 0
 weight-constr2-t0-vN10: - 524288 Y[N10] + Z[N10] <= 0
 weight-constr3-t0-vN10: Y[N10] - Z[N10] <= 0
 weight-constr2-t0-vN13: - 524288 Y[N13] + Z[N13] <= 0
 weight-constr3-t0-vN13: Y[N13] - Z[N13] <= 0
 weight-constr2-t0-vN1: - 524288 Y[N1] + Z[N1] <= 0
 weight-constr3-t0-vN1: Y[N1] - Z[N1] <= 0
 weight-constr2-t1-vN11: - 524288 Y[N11] + Z[N11] <= 0
 weight-constr3-t1-vN11: Y[N11] - Z[N11] <= 0
 weight-constr2-t1-vN7: - 524288 Y[N7] + Z[N7] <= 0
 weight-constr3-t1-vN7: Y[N7] - Z[N7] <= 0
 weight-constr2-t1-vN8: - 524288 Y[N8] + Z[N8] <= 0
 weight-constr3-t1-vN8: Y[N8] - Z[N8] <= 0
 weight-constr2-t1-vN4: - 524288 Y[N4] + Z[N4] <= 0
 weight-constr3-t1-vN4: Y[N4] - Z[N4] <= 0
 weight-constr2-t1-vN2: - 524288 Y[N2] + Z[N2] <= 0
 weight-constr3-t1-vN2: Y[N2] - Z[N2] <= 0
 weight-constr2-t1-vN15: - 524288 Y[N15] + Z[N15] <= 0
 weight-constr3-t1-vN15: Y[N15] - Z[N15] <= 0
 weight-constr2-t1-vN5: - 524288 Y[N5] + Z[N5] <= 0
 weight-constr3-t1-vN5: Y[N5] - Z[N5] <= 0
 weight-constr2-t1-vN6: - 524288 Y[N6] + Z[N6] <= 0
 weight-constr3-t1-vN6: Y[N6] - Z[N6] <= 0
 weight-constr2-t1-vN3: - 524288 Y[N3] + Z[N3] <= 0
 weight-constr3-t1-vN3: Y[N3] - Z[N3] <= 0
 weight-constr2-t1-vN12: - 524288 Y[N12] + Z[N12] <= 0
 weight-constr3-t1-vN12: Y[N12] - Z[N12] <= 0
 weight-constr2-t1-vN14: - 524288 Y[N14] + Z[N14] <= 0
 weight-constr3-t1-vN14: Y[N14] - Z[N14] <= 0
 weight-constr2-t1-vN9: - 524288 Y[N9] + Z[N9] <= 0
 weight-constr3-t1-vN9: Y[N9] - Z[N9] <= 0
 weight-constr2-t1-vN10: - 524288 Y[N10] + Z[N10] <= 0
 weight-constr3-t1-vN10: Y[N10] - Z[N10] <= 0
 weight-constr2-t1-vN13: - 524288 Y[N13] + Z[N13] <= 0
 weight-constr3-t1-vN13: Y[N13] - Z[N13] <= 0
 weight-constr2-t1-vN1: - 524288 Y[N1] + Z[N1] <= 0
 weight-constr3-t1-vN1: Y[N1] - Z[N1] <= 0
 weight-constr2-t2-vN11: - 524288 Y[N11] + Z[N11] <= 0
 weight-constr3-t2-vN11: Y[N11] - Z[N11] <= 0
 weight-constr2-t2-vN7: - 524288 Y[N7] + Z[N7] <= 0
 weight-constr3-t2-vN7: Y[N7] - Z[N7] <= 0
 weight-constr2-t2-vN8: - 524288 Y[N8] + Z[N8] <= 0
 weight-constr3-t2-vN8: Y[N8] - Z[N8] <= 0
 weight-constr2-t2-vN4: - 524288 Y[N4] + Z[N4] <= 0
 weight-constr3-t2-vN4: Y[N4] - Z[N4] <= 0
 weight-constr2-t2-vN2: - 524288 Y[N2] + Z[N2] <= 0
 weight-constr3-t2-vN2: Y[N2] - Z[N2] <= 0
 weight-constr2-t2-vN15: - 524288 Y[N15] + Z[N15] <= 0
 weight-constr3-t2-vN15: Y[N15] - Z[N15] <= 0
 weight-constr2-t2-vN5: - 524288 Y[N5] + Z[N5] <= 0
 weight-constr3-t2-vN5: Y[N5] - Z[N5] <= 0
 weight-constr2-t2-vN6: - 524288 Y[N6] + Z[N6] <= 0
 weight-constr3-t2-vN6: Y[N6] - Z[N6] <= 0
 weight-constr2-t2-vN3: - 524288 Y[N3] + Z[N3] <= 0
 weight-constr3-t2-vN3: Y[N3] - Z[N3] <= 0
 weight-constr2-t2-vN12: - 524288 Y[N12] + Z[N12] <= 0
 weight-constr3-t2-vN12: Y[N12] - Z[N12] <= 0
 weight-constr2-t2-vN14: - 524288 Y[N14] + Z[N14] <= 0
 weight-constr3-t2-vN14: Y[N14] - Z[N14] <= 0
 weight-constr2-t2-vN9: - 524288 Y[N9] + Z[N9] <= 0
 weight-constr3-t2-vN9: Y[N9] - Z[N9] <= 0
 weight-constr2-t2-vN10: - 524288 Y[N10] + Z[N10] <= 0
 weight-constr3-t2-vN10: Y[N10] - Z[N10] <= 0
 weight-constr2-t2-vN13: - 524288 Y[N13] + Z[N13] <= 0
 weight-constr3-t2-vN13: Y[N13] - Z[N13] <= 0
 weight-constr2-t2-vN1: - 524288 Y[N1] + Z[N1] <= 0
 weight-constr3-t2-vN1: Y[N1] - Z[N1] <= 0
 weight-constr2-t3-vN11: - 524288 Y[N11] + Z[N11] <= 0
 weight-constr3-t3-vN11: Y[N11] - Z[N11] <= 0
 weight-constr2-t3-vN7: - 524288 Y[N7] + Z[N7] <= 0
 weight-constr3-t3-vN7: Y[N7] - Z[N7] <= 0
 weight-constr2-t3-vN8: - 524288 Y[N8] + Z[N8] <= 0
 weight-constr3-t3-vN8: Y[N8] - Z[N8] <= 0
 weight-constr2-t3-vN4: - 524288 Y[N4] + Z[N4] <= 0
 weight-constr3-t3-vN4: Y[N4] - Z[N4] <= 0
 weight-constr2-t3-vN2: - 524288 Y[N2] + Z[N2] <= 0
 weight-constr3-t3-vN2: Y[N2] - Z[N2] <= 0
 weight-constr2-t3-vN15: - 524288 Y[N15] + Z[N15] <= 0
 weight-constr3-t3-vN15: Y[N15] - Z[N15] <= 0
 weight-constr2-t3-vN5: - 524288 Y[N5] + Z[N5] <= 0
 weight-constr3-t3-vN5: Y[N5] - Z[N5] <= 0
 weight-constr2-t3-vN6: - 524288 Y[N6] + Z[N6] <= 0
 weight-constr3-t3-vN6: Y[N6] - Z[N6] <= 0
 weight-constr2-t3-vN3: - 524288 Y[N3] + Z[N3] <= 0
 weight-constr3-t3-vN3: Y[N3] - Z[N3] <= 0
 weight-constr2-t3-vN12: - 524288 Y[N12] + Z[N12] <= 0
 weight-constr3-t3-vN12: Y[N12] - Z[N12] <= 0
 weight-constr2-t3-vN14: - 524288 Y[N14] + Z[N14] <= 0
 weight-constr3-t3-vN14: Y[N14] - Z[N14] <= 0
 weight-constr2-t3-vN9: - 524288 Y[N9] + Z[N9] <= 0
 weight-constr3-t3-vN9: Y[N9] - Z[N9] <= 0
 weight-constr2-t3-vN10: - 524288 Y[N10] + Z[N10] <= 0
 weight-constr3-t3-vN10: Y[N10] - Z[N10] <= 0
 weight-constr2-t3-vN13: - 524288 Y[N13] + Z[N13] <= 0
 weight-constr3-t3-vN13: Y[N13] - Z[N13] <= 0
 weight-constr2-t3-vN1: - 524288 Y[N1] + Z[N1] <= 0
 weight-constr3-t3-vN1: Y[N1] - Z[N1] <= 0
 weight-constr2-t4-vN11: - 524288 Y[N11] + Z[N11] <= 0
 weight-constr3-t4-vN11: Y[N11] - Z[N11] <= 0
 weight-constr2-t4-vN7: - 524288 Y[N7] + Z[N7] <= 0
 weight-constr3-t4-vN7: Y[N7] - Z[N7] <= 0
 weight-constr2-t4-vN8: - 524288 Y[N8] + Z[N8] <= 0
 weight-constr3-t4-vN8: Y[N8] - Z[N8] <= 0
 weight-constr2-t4-vN4: - 524288 Y[N4] + Z[N4] <= 0
 weight-constr3-t4-vN4: Y[N4] - Z[N4] <= 0
 weight-constr2-t4-vN2: - 524288 Y[N2] + Z[N2] <= 0
 weight-constr3-t4-vN2: Y[N2] - Z[N2] <= 0
 weight-constr2-t4-vN15: - 524288 Y[N15] + Z[N15] <= 0
 weight-constr3-t4-vN15: Y[N15] - Z[N15] <= 0
 weight-constr2-t4-vN5: - 524288 Y[N5] + Z[N5] <= 0
 weight-constr3-t4-vN5: Y[N5] - Z[N5] <= 0
 weight-constr2-t4-vN6: - 524288 Y[N6] + Z[N6] <= 0
 weight-constr3-t4-vN6: Y[N6] - Z[N6] <= 0
 weight-constr2-t4-vN3: - 524288 Y[N3] + Z[N3] <= 0
 weight-constr3-t4-vN3: Y[N3] - Z[N3] <= 0
 weight-constr2-t4-vN12: - 524288 Y[N12] + Z[N12] <= 0
 weight-constr3-t4-vN12: Y[N12] - Z[N12] <= 0
 weight-constr2-t4-vN14: - 524288 Y[N14] + Z[N14] <= 0
 weight-constr3-t4-vN14: Y[N14] - Z[N14] <= 0
 weight-constr2-t4-vN9: - 524288 Y[N9] + Z[N9] <= 0
 weight-constr3-t4-vN9: Y[N9] - Z[N9] <= 0
 weight-constr2-t4-vN10: - 524288 Y[N10] + Z[N10] <= 0
 weight-constr3-t4-vN10: Y[N10] - Z[N10] <= 0
 weight-constr2-t4-vN13: - 524288 Y[N13] + Z[N13] <= 0
 weight-constr3-t4-vN13: Y[N13] - Z[N13] <= 0
 weight-constr2-t4-vN1: - 524288 Y[N1] + Z[N1] <= 0
 weight-constr3-t4-vN1: Y[N1] - Z[N1] <= 0
 weight-constr2-t5-vN11: - 524288 Y[N11] + Z[N11] <= 0
 weight-constr3-t5-vN11: Y[N11] - Z[N11] <= 0
 weight-constr2-t5-vN7: - 524288 Y[N7] + Z[N7] <= 0
 weight-constr3-t5-vN7: Y[N7] - Z[N7] <= 0
 weight-constr2-t5-vN8: - 524288 Y[N8] + Z[N8] <= 0
 weight-constr3-t5-vN8: Y[N8] - Z[N8] <= 0
 weight-constr2-t5-vN4: - 524288 Y[N4] + Z[N4] <= 0
 weight-constr3-t5-vN4: Y[N4] - Z[N4] <= 0
 weight-constr2-t5-vN2: - 524288 Y[N2] + Z[N2] <= 0
 weight-constr3-t5-vN2: Y[N2] - Z[N2] <= 0
 weight-constr2-t5-vN15: - 524288 Y[N15] + Z[N15] <= 0
 weight-constr3-t5-vN15: Y[N15] - Z[N15] <= 0
 weight-constr2-t5-vN5: - 524288 Y[N5] + Z[N5] <= 0
 weight-constr3-t5-vN5: Y[N5] - Z[N5] <= 0
 weight-constr2-t5-vN6: - 524288 Y[N6] + Z[N6] <= 0
 weight-constr3-t5-vN6: Y[N6] - Z[N6] <= 0
 weight-constr2-t5-vN3: - 524288 Y[N3] + Z[N3] <= 0
 weight-constr3-t5-vN3: Y[N3] - Z[N3] <= 0
 weight-constr2-t5-vN12: - 524288 Y[N12] + Z[N12] <= 0
 weight-constr3-t5-vN12: Y[N12] - Z[N12] <= 0
 weight-constr2-t5-vN14: - 524288 Y[N14] + Z[N14] <= 0
 weight-constr3-t5-vN14: Y[N14] - Z[N14] <= 0
 weight-constr2-t5-vN9: - 524288 Y[N9] + Z[N9] <= 0
 weight-constr3-t5-vN9: Y[N9] - Z[N9] <= 0
 weight-constr2-t5-vN10: - 524288 Y[N10] + Z[N10] <= 0
 weight-constr3-t5-vN10: Y[N10] - Z[N10] <= 0
 weight-constr2-t5-vN13: - 524288 Y[N13] + Z[N13] <= 0
 weight-constr3-t5-vN13: Y[N13] - Z[N13] <= 0
 weight-constr2-t5-vN1: - 524288 Y[N1] + Z[N1] <= 0
 weight-constr3-t5-vN1: Y[N1] - Z[N1] <= 0
Bounds
 Z[N11] <= 524288
 Z[N7] <= 524288
 Z[N8] <= 524288
 Z[N4] <= 524288
 Z[N2] <= 524288
 Z[N15] <= 524288
 Z[N5] <= 524288
 Z[N6] <= 524288
 Z[N3] <= 524288
 Z[N12] <= 524288
 Z[N14] <= 524288
 Z[N9] <= 524288
 Z[N10] <= 524288
 Z[N13] <= 524288
 Z[N1] <= 524288
 Z[N11] <= 524288
 Z[N7] <= 524288
 Z[N8] <= 524288
 Z[N4] <= 524288
 Z[N2] <= 524288
 Z[N15] <= 524288
 Z[N5] <= 524288
 Z[N6] <= 524288
 Z[N3] <= 524288
 Z[N12] <= 524288
 Z[N14] <= 524288
 Z[N9] <= 524288
 Z[N10] <= 524288
 Z[N13] <= 524288
 Z[N1] <= 524288
 Z[N11] <= 524288
 Z[N7] <= 524288
 Z[N8] <= 524288
 Z[N4] <= 524288
 Z[N2] <= 524288
 Z[N15] <= 524288
 Z[N5] <= 524288
 Z[N6] <= 524288
 Z[N3] <= 524288
 Z[N12] <= 524288
 Z[N14] <= 524288
 Z[N9] <= 524288
 Z[N10] <= 524288
 Z[N13] <= 524288
 Z[N1] <= 524288
 Z[N11] <= 524288
 Z[N7] <= 524288
 Z[N8] <= 524288
 Z[N4] <= 524288
 Z[N2] <= 524288
 Z[N15] <= 524288
 Z[N5] <= 524288
 Z[N6] <= 524288
 Z[N3] <= 524288
 Z[N12] <= 524288
 Z[N14] <= 524288
 Z[N9] <= 524288
 Z[N10] <= 524288
 Z[N13] <= 524288
 Z[N1] <= 524288
 Z[N11] <= 524288
 Z[N7] <= 524288
 Z[N8] <= 524288
 Z[N4] <= 524288
 Z[N2] <= 524288
 Z[N15] <= 524288
 Z[N5] <= 524288
 Z[N6] <= 524288
 Z[N3] <= 524288
 Z[N12] <= 524288
 Z[N14] <= 524288
 Z[N9] <= 524288
 Z[N10] <= 524288
 Z[N13] <= 524288
 Z[N1] <= 524288
 Z[N11] <= 524288
 Z[N7] <= 524288
 Z[N8] <= 524288
 Z[N4] <= 524288
 Z[N2] <= 524288
 Z[N15] <= 524288
 Z[N5] <= 524288
 Z[N6] <= 524288
 Z[N3] <= 524288
 Z[N12] <= 524288
 Z[N14] <= 524288
 Z[N9] <= 524288
 Z[N10] <= 524288
 Z[N13] <= 524288
 Z[N1] <= 524288
Binaries
 X[N6,N1] X[N1,N6] X[N7,N1] X[N1,N7] X[N8,N1] X[N1,N8] X[N3,N2] X[N2,N3]
 X[N5,N2] X[N2,N5] X[N6,N2] X[N2,N6] X[N5,N3] X[N3,N5] X[N8,N3] X[N3,N8]
 X[N5,N4] X[N4,N5] X[N6,N4] X[N4,N6] X[N13,N6] X[N6,N13] X[N10,N7]
 X[N7,N10] X[N14,N7] X[N7,N14] X[N9,N8] X[N8,N9] X[N15,N8] X[N8,N15]
 X[N10,N9] X[N9,N10] X[N12,N9] X[N9,N12] X[N15,N9] X[N9,N15] X[N12,N10]
 X[N10,N12] X[N13,N11] X[N11,N13] X[N14,N11] X[N11,N14] X[N14,N13]
 X[N13,N14] Y[N11] Y[N7] Y[N8] Y[N4] Y[N2] Y[N15] Y[N5] Y[N6] Y[N3] Y[N12]
 Y[N14] Y[N9] Y[N10] Y[N13] Y[N1] X[N6,N1] X[N1,N6] X[N7,N1] X[N1,N7]
 X[N8,N1] X[N1,N8] X[N3,N2] X[N2,N3] X[N5,N2] X[N2,N5] X[N6,N2] X[N2,N6]
 X[N5,N3] X[N3,N5] X[N8,N3] X[N3,N8] X[N5,N4] X[N4,N5] X[N6,N4] X[N4,N6]
 X[N13,N6] X[N6,N13] X[N10,N7] X[N7,N10] X[N14,N7] X[N7,N14] X[N9,N8]
 X[N8,N9] X[N15,N8] X[N8,N15] X[N10,N9] X[N9,N10] X[N12,N9] X[N9,N12]
 X[N15,N9] X[N9,N15] X[N12,N10] X[N10,N12] X[N13,N11] X[N11,N13] X[N14,N11]
 X[N11,N14] X[N14,N13] X[N13,N14] Y[N11] Y[N7] Y[N8] Y[N4] Y[N2] Y[N15]
 Y[N5] Y[N6] Y[N3] Y[N12] Y[N14] Y[N9] Y[N10] Y[N13] Y[N1] X[N6,N1]
 X[N1,N6] X[N7,N1] X[N1,N7] X[N8,N1] X[N1,N8] X[N3,N2] X[N2,N3] X[N5,N2]
 X[N2,N5] X[N6,N2] X[N2,N6] X[N5,N3] X[N3,N5] X[N8,N3] X[N3,N8] X[N5,N4]
 X[N4,N5] X[N6,N4] X[N4,N6] X[N13,N6] X[N6,N13] X[N10,N7] X[N7,N10]
 X[N14,N7] X[N7,N14] X[N9,N8] X[N8,N9] X[N15,N8] X[N8,N15] X[N10,N9]
 X[N9,N10] X[N12,N9] X[N9,N12] X[N15,N9] X[N9,N15] X[N12,N10] X[N10,N12]
 X[N13,N11] X[N11,N13] X[N14,N11] X[N11,N14] X[N14,N13] X[N13,N14] Y[N11]
 Y[N7] Y[N8] Y[N4] Y[N2] Y[N15] Y[N5] Y[N6] Y[N3] Y[N12] Y[N14] Y[N9]
 Y[N10] Y[N13] Y[N1] X[N6,N1] X[N1,N6] X[N7,N1] X[N1,N7] X[N8,N1] X[N1,N8]
 X[N3,N2] X[N2,N3] X[N5,N2] X[N2,N5] X[N6,N2] X[N2,N6] X[N5,N3] X[N3,N5]
 X[N8,N3] X[N3,N8] X[N5,N4] X[N4,N5] X[N6,N4] X[N4,N6] X[N13,N6] X[N6,N13]
 X[N10,N7] X[N7,N10] X[N14,N7] X[N7,N14] X[N9,N8] X[N8,N9] X[N15,N8]
 X[N8,N15] X[N10,N9] X[N9,N10] X[N12,N9] X[N9,N12] X[N15,N9] X[N9,N15]
 X[N12,N10] X[N10,N12] X[N13,N11] X[N11,N13] X[N14,N11] X[N11,N14]
 X[N14,N13] X[N13,N14] Y[N11] Y[N7] Y[N8] Y[N4] Y[N2] Y[N15] Y[N5] Y[N6]
 Y[N3] Y[N12] Y[N14] Y[N9] Y[N10] Y[N13] Y[N1] X[N6,N1] X[N1,N6] X[N7,N1]
 X[N1,N7] X[N8,N1] X[N1,N8] X[N3,N2] X[N2,N3] X[N5,N2] X[N2,N5] X[N6,N2]
 X[N2,N6] X[N5,N3] X[N3,N5] X[N8,N3] X[N3,N8] X[N5,N4] X[N4,N5] X[N6,N4]
 X[N4,N6] X[N13,N6] X[N6,N13] X[N10,N7] X[N7,N10] X[N14,N7] X[N7,N14]
 X[N9,N8] X[N8,N9] X[N15,N8] X[N8,N15] X[N10,N9] X[N9,N10] X[N12,N9]
 X[N9,N12] X[N15,N9] X[N9,N15] X[N12,N10] X[N10,N12] X[N13,N11] X[N11,N13]
 X[N14,N11] X[N11,N14] X[N14,N13] X[N13,N14] Y[N11] Y[N7] Y[N8] Y[N4] Y[N2]
 Y[N15] Y[N5] Y[N6] Y[N3] Y[N12] Y[N14] Y[N9] Y[N10] Y[N13] Y[N1] X[N6,N1]
 X[N1,N6] X[N7,N1] X[N1,N7] X[N8,N1] X[N1,N8] X[N3,N2] X[N2,N3] X[N5,N2]
 X[N2,N5] X[N6,N2] X[N2,N6] X[N5,N3] X[N3,N5] X[N8,N3] X[N3,N8] X[N5,N4]
 X[N4,N5] X[N6,N4] X[N4,N6] X[N13,N6] X[N6,N13] X[N10,N7] X[N7,N10]
 X[N14,N7] X[N7,N14] X[N9,N8] X[N8,N9] X[N15,N8] X[N8,N15] X[N10,N9]
 X[N9,N10] X[N12,N9] X[N9,N12] X[N15,N9] X[N9,N15] X[N12,N10] X[N10,N12]
 X[N13,N11] X[N11,N13] X[N14,N11] X[N11,N14] X[N14,N13] X[N13,N14] Y[N11]
 Y[N7] Y[N8] Y[N4] Y[N2] Y[N15] Y[N5] Y[N6] Y[N3] Y[N12] Y[N14] Y[N9]
 Y[N10] Y[N13] Y[N1]
End
