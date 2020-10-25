import seqfhfit as seq
import numpy as np
import gvar as gv

# 7 states for 2pt and 3pt, 4 states for sum and fh,
# t = [4, 10] for all
X = np.arange(4, 11)

# results from jin chen's code
# for 3pt this is the time ordering:
# 'pt3_A3': tsep: [array([ 4,  4,  5,  5,  6,  6,  6,  7,  7,  7,  8,  8,  8,  8,  9,  9,  9, 9, 10, 10, 10, 10, 10]),
#           tau: array([1, 2, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 5])]
result = {'pt2': np.array(
    [4.68440257e-08, 2.29203824e-08, 1.14692352e-08, 5.81066551e-09, 2.96513323e-09, 1.51975678e-09, 7.81146394e-10]),
    'pt3_A3': np.array(
        [5.83510860e-08, 5.28737617e-08, 2.89829869e-08, 2.60435153e-08, 1.46751859e-08, 1.32056589e-08, 1.30813079e-08, 7.50191550e-09, 6.76977488e-09, 6.70768143e-09, 3.85431159e-09, 3.48707511e-09, 3.46166678e-09, 3.46167078e-09, 1.98583328e-09, 1.80019141e-09, 1.79043923e-09, 1.79341886e-09, 1.02486321e-09, 9.30410326e-10, 9.26743015e-10, 9.29848275e-10, 9.31333486e-10]),
    'pt3_V4': np.array(
        [5.29564757e-08, 4.93608640e-08, 2.58553820e-08, 2.38479022e-08, 1.29270650e-08, 1.18849768e-08, 1.17635678e-08, 6.54612196e-09, 6.01308962e-09, 5.93287695e-09, 3.33951012e-09, 3.06709078e-09, 3.02352814e-09, 3.01384490e-09, 1.71149653e-09, 1.57189082e-09, 1.54937260e-09, 1.54295815e-09, 8.79770741e-10, 8.08027035e-10, 7.96521788e-10, 7.93084149e-10, 7.92312131e-10]),
    'sum_A3': np.array([1.16476226, 1.19240865, 1.21439542, 1.231113, 1.24352231, 1.2526004,  1.2591593]),
    'sum_V4': np.array([1.00699138, 1.00767945, 1.00876142, 1.01038971, 1.01241361, 1.01462981, 1.0168711]),
    'fh_ss_A3': np.array([1.24856479, 1.24433982, 1.24548882, 1.2492698, 1.25388609, 1.25836057, 1.26224102]),
    'fh_ss_V4': np.array([0.9987932, 1.00885136, 1.01421799, 1.01741907, 1.01966673, 1.02148297, 1.02307095])}

init = {"nt": 50, "strategy": "does not matter"}
ff = seq.FitFunction(init)

# for 2pt and 3pt
p = {'e0': 0.6569980804940844, 'log(e1)': -1.7925646342559884, 'z0': 0.0007322962713869645, 'z1': 0.0003192673971732764, 'z0_ps': 0.0027599222530335322, 'z1_ps': 0.0014311182913749262, 'log(e2)': -0.9846627109449616, 'z2': 0.0005867499848728282, 'z2_ps': 0.004925016671055348, 'log(e3)': -1.433263239979172, 'z3': 0.0005857677045485657, 'z3_ps': 0.003743285963156695, 'log(e4)': -0.8783809400029784, 'z4': 0.0007358919366326382, 'z4_ps': 0.004101609401627379,  'log(e5)': -0.5990846646746402, 'z5': 0.0007209782743180652, 'z5_ps': 0.002475801393447248, 'log(e6)': 0.22087872995626118, 'z6': 0.001068604493125576, 'z6_ps': 0.0015274763367897285, 'a00': 1.2605222674623597, 'v00': 1.034619576802522, 'a10': -0.019928993130168368, 'v10': -0.08257541644904179, 'a11': 0.2919843386539549, 'v11': 1.0511358174118686, 'a20': -0.4390643954205519, 'v20': 0.12174272715835381, 'a21': 0.3620476349009373, 'v21': 0.1833710353298747, 'a22': 0.6361286846700358, 'v22': 0.9421169783172731, 'a30': 0.34375532136084036, 'v30': -0.2905312734758731, 'a31': 0.014756209659542678, 'v31': -0.08191836353816004, 'a32': 0.425140113428051, 'v32': -0.25852654276212567, 'a33': 0.0761691021816558, 'v33': 0.9768394026150242, 'a40': 0.3228590485478712, 'v40': 0.6600701077522847, 'a41': -0.20601569608454806,'v41': 0.120214900669482, 'a42': 0.1357338625113631, 'v42': 0.9759940383034039, 'a43': -0.008233363922618956, 'v43': -0.589533998863517, 'a44': -0.08441553733853056, 'v44': 1.0169700028968167, 'a50': -0.2895742440269932, 'v50': -0.6859614174750054, 'a51': 0.007091833315171467, 'v51': -0.4825540980298498, 'a52': 0.24375784517752733,'v52': -0.17342902031227236, 'a53': 0.3622536648264607, 'v53': -0.8199613170807609, 'a54': 0.2883915685483757, 'v54': 0.653503183904798, 'a55': 0.1388100787683668, 'v55': 1.0116094504724715, 'a60': 1.5571138917093639, 'v60': 1.2648036160630534, 'a61': 0.6498508782851965, 'v61': 0.6630205376714541, 'a62': 0.3760553976876156, 'v62': 0.7303543467645489, 'a63': 0.4681752983016586, 'v63': 0.22351478241509143, 'a64': 0.43636684770879647, 'v64': 0.7836973303190226, 'a65': 0.1989987756810002, 'v65': 0.3843114716519219, 'a66': 0.014721329494868908, 'v66': 0.9997916964382514}
p = gv.BufferDict(p)

ff.strategy = dict()
ff.strategy["nstates"] = 7

# check 2pt
twopt = ff.twopoint(X, p)
diff_twopt = (result["pt2"] - twopt) / result["pt2"]
print("two point relative difference")
print(diff_twopt)

# check 3pt
# the ceiling function cuts the data in half, which assumes the data is folded vs. tau
from math import ceil

a3 = []
v4 = []
for Xi in X:
    t = np.arange(1, Xi)
    ff.corr = [False, "a", Xi]
    r = ff.threepoint(t, p)
    a3.extend(r[:ceil(len(r) / 2)])
    ff.corr = [False, "v", Xi]
    r = ff.threepoint(t, p)
    v4.extend(r[:ceil(len(r) / 2)])
print("three point relative different")
print("A3")
print((result["pt3_A3"] - a3) / result["pt3_A3"])
print("V4")
print((result["pt3_V4"] - v4) / result["pt3_V4"])

# for sum and fh
# take log(e_fh) -> log(e3) and z_fh_ss -> z3
# fh_a0,1,2,3 -> a30, a31, a32, a33
# fh_v0,1,2,3 -> v30, v31, v32, v33
# dn_ss_A3 -> cn
# dn_ss_V4 -> dn
p = {'e0': 0.6569980804940844, 'log(e1)': -1.7925646342559884, 'z0': 0.0007322962713869645, 'z1': 0.0003192673971732764, 'z0_ps': 0.0027599222530335322, 'z1_ps': 0.0014311182913749262, 'log(e2)': -0.9846627109449616, 'z2': 0.0005867499848728282, 'z2_ps': 0.004925016671055348,  'log(e3)': -1.3401326343142712, 'z3': 0.0005814368336925668, 'z3_ps': 0.002254903914697272, 'a00': 1.2605222674623597, 'v00': 1.034619576802522, 'a10': -0.019928993130168368, 'v10': -0.08257541644904179, 'a11': 0.2919843386539549, 'v11': 1.0511358174118686, 'a20': -0.4390643954205519, 'v20': 0.12174272715835381, 'a21': 0.3620476349009373, 'v21': 0.1833710353298747, 'a22': 0.6361286846700358, 'v22': 0.9421169783172731,
'a30': 0.4071791650540326, 'v30': -0.1988939247131013, 'a31': -0.11848922973337624, 'v31': -0.162396828097687, 'a32': 0.4755531742816833, 'v32': 0.07142936237560996, 'a33': 0.31855072576988946, 'v33': 1.0119119986917513, 'c0': -1.9243967556021374e-06, 'd0': 9.802455759543743e-07, 'c1': -3.4653052773960966e-07, 'd1': 5.270038139207717e-08, 'c2': -2.43881113104576e-06, 'd2': 7.455135470437805e-07, 'c3': -1.2715045483291115e-06, 'd3': 9.716182722893693e-07}
p = gv.BufferDict(p)

# check sum
ff.strategy["nstates"] = 4
ff.strategy["fstates"] = 4
ff.corr = [False, "sma"]
smA = ff.fh(X, p, "sm")
ff.corr = [False, "smv"]
smV = ff.fh(X, p, "sm")

print("sum relative diff")
print("A3")
print((result["sum_A3"] - smA) / result["sum_A3"])
print("V4")
print((result["sum_V4"] - smV) / result["sum_V4"])

# check fh
ff.corr = [False, "fha"]
fhA = ff.fh(X, p, "fh")
ff.corr = [False, "fhv"]
fhV = ff.fh(X, p, "fh")

print("fh relative diff")
print("A3")
print((result["fh_ss_A3"] - fhA) / result["fh_ss_A3"])
print("V4")
print((result["fh_ss_V4"] - fhV) / result["fh_ss_V4"])