def read_parameters(ensemble):
    parameters = dict()

    # metadata
    parameters["ensemble"] = ensemble
    if ensemble == "a09m310":
        parameters["h5filename"] = "a09m310_e_avg.h5"
        parameters["h5prepath"] = "gf1p0_w3p5_n45_M51p1_L56_a1p5"
        parameters["mval"] = 0.00951
        parameters["tsep"] = range(2, 15)
    elif ensemble == "a12m310":
        parameters["h5filename"] = "a12m310_a_avg.h5"
        parameters["h5prepath"] = "gf1p0_w3p0_n30_M51p2_L58_a1p5"
        parameters["mval"] = 0.0126
        parameters["tsep"] = range(4, 12)
    elif ensemble == "a15m310":
        parameters["h5filename"] = "a13m310L_a_avg.h5"
        parameters["h5prepath"] = "gf1p0_w3p0_n30_M51p3_L512_a1p5"
        parameters["mval"] = 0.0158
        parameters["tsep"] = range(2, 11)
    elif ensemble == "a09m135":
        parameters["h5filename"] = "a09m135_s_avg_srcs0-15.h5"
        parameters["h5prepath"] = "gf1p0_w3p5_n45_M51p1_L512_a2p0"
        parameters["mval"] = 0.00152
        parameters["tsep"] = range(3, 13)
    elif ensemble == "a12m135":
        parameters["h5filename"] = "a12m130_a_avg_srcs0-31.h5"
        parameters["h5prepath"] = "gf1p0_w3p0_n30_M51p2_L520_a3p0"
        parameters["mval"] = 0.00195
        parameters["tsep"] = range(3, 12)

    # fit strategy
    # nstates defines number of states for all correlators except FH
    # fstates defines number of FH states
    # if nstates != fstates, the code will spawn a new trashcan state for FH w/ same (uncorr) prior
    # nstate >= fstate otherwise the priors will be set up incorrectly. This is the option that makes sense anyways.
    # sequence is a list of lists defining the set of chained simultaneous fits
    # sequence defines tmin and tmax for seq3pt as well
    parameters["strategy"] = dict()
    p = parameters["strategy"]
    if ensemble == "a09m310":
        p["nstates"] = 6
        p["fstates"] = 2
        p["sequence"] = [
            [
                ("proton", "twopt"),
                ("proton", "A3", 2),
                ("proton", "A3", 3),
                ("proton", "A3", 4),
                ("proton", "A3", 5),
                ("proton", "A3", 6),
                ("proton", "A3", 7),
                ("proton", "A3", 8),
                ("proton", "A3", 9),
                ("proton", "A3", 10),
                ("proton", "V4", 2),
                ("proton", "V4", 3),
                ("proton", "V4", 4),
                ("proton", "V4", 5),
                ("proton", "V4", 6),
                ("proton", "V4", 7),
                ("proton", "V4", 8),
                ("proton", "V4", 9),
                ("proton", "V4", 10),
                ("proton", "smA3"),
                ("proton", "smV4")
            ]
        ]
        #p["sequence"] = [[("proton", "twopt"), ("proton", "fhA3"), ("proton", "fhV4")]]
        #p["sequence"] = [[("proton", "twopt")]]
        # trange for seq3pt is for insertion-sink separation time
        p["trange"] = dict()
        t = p["trange"]
        t[("proton", "twopt")] = [4, 15]
        t[("proton", "smA3")] = [4, 8]
        t[("proton", "smV4")] = [4, 8]
        t[("proton", "A3", 2)] = [1, 1]
        t[("proton", "A3", 3)] = [1, 2]
        t[("proton", "A3", 4)] = [1, 3]
        t[("proton", "A3", 5)] = [1, 4]
        t[("proton", "A3", 6)] = [1, 5]
        t[("proton", "A3", 7)] = [1, 6]
        t[("proton", "A3", 8)] = [1, 7]
        t[("proton", "A3", 9)] = [1, 8]
        t[("proton", "A3", 10)] = [1, 9]
        t[("proton", "V4", 2)] = [1, 1]
        t[("proton", "V4", 3)] = [1, 2]
        t[("proton", "V4", 4)] = [1, 3]
        t[("proton", "V4", 5)] = [1, 4]
        t[("proton", "V4", 6)] = [1, 5]
        t[("proton", "V4", 7)] = [1, 6]
        t[("proton", "V4", 8)] = [1, 7]
        t[("proton", "V4", 9)] = [1, 8]
        t[("proton", "V4", 10)] = [1, 9]
    elif ensemble == "a12m310":
        pass
    elif ensemble == "a15m310":
        pass
    elif ensemble == "a09m135":
        p["nstates"] = 8
        p["fstates"] = 4
        p["sequence"] = [
            [
                ("proton", "twopt"),
                #("proton", "A3", 3),
                ("proton", "A3", 4),
                ("proton", "A3", 5),
                ("proton", "A3", 6),
                ("proton", "A3", 7),
                ("proton", "A3", 8),
                ("proton", "A3", 9),
                ("proton", "A3", 10),
                ("proton", "A3", 11),
                ("proton", "A3", 12),
                #("proton", "V4", 3),
                ("proton", "V4", 4),
                ("proton", "V4", 5),
                ("proton", "V4", 6),
                ("proton", "V4", 7),
                ("proton", "V4", 8),
                ("proton", "V4", 9),
                ("proton", "V4", 10),
                ("proton", "V4", 11),
                ("proton", "V4", 12),
                ("proton", "fhA3"),
                ("proton", "fhV4")
            ]
        ]
        #p["sequence"] = [[("proton", "twopt"), ("proton", "fhA3"), ("proton", "fhV4")]]
        #p["sequence"] = [[("proton", "twopt")]]
        # trange for seq3pt is for insertion-sink separation time
        p["trange"] = dict()
        t = p["trange"]
        t[("proton", "twopt")] = [4, 15]
        t[("proton", "fhA3")] = [4, 11]
        t[("proton", "fhV4")] = [4, 11]
        t[("proton", "A3", 3)] = [1, 2]
        t[("proton", "A3", 4)] = [1, 3]
        t[("proton", "A3", 5)] = [1, 4]
        t[("proton", "A3", 6)] = [1, 5]
        t[("proton", "A3", 7)] = [1, 6]
        t[("proton", "A3", 8)] = [1, 7]
        t[("proton", "A3", 9)] = [1, 8]
        t[("proton", "A3", 10)] = [1, 9]
        t[("proton", "A3", 11)] = [1, 10]
        t[("proton", "A3", 12)] = [1, 11]
        t[("proton", "V4", 3)] = [1, 2]
        t[("proton", "V4", 4)] = [1, 3]
        t[("proton", "V4", 5)] = [1, 4]
        t[("proton", "V4", 6)] = [1, 5]
        t[("proton", "V4", 7)] = [1, 6]
        t[("proton", "V4", 8)] = [1, 7]
        t[("proton", "V4", 9)] = [1, 8]
        t[("proton", "V4", 10)] = [1, 9]
        t[("proton", "V4", 11)] = [1, 10]
        t[("proton", "V4", 12)] = [1, 11]
    elif ensemble == "a12m135":
        p["nstates"] = 4
        p["fstates"] = 3
        p["sequence"] = [
            [
                ("proton", "twopt"),
                ("proton", "smA3"),
                ("proton", "smV4"),
                ("proton", "fhA3"),
                ("proton", "fhV4"),
                #("proton", "A3", 3),
                ("proton", "A3", 4),
                ("proton", "A3", 5),
                ("proton", "A3", 6),
                ("proton", "A3", 7),
                ("proton", "A3", 8),
                ("proton", "A3", 9),
                ("proton", "A3", 10),
                #("proton", "A3", 11),
                #("proton", "V4", 3),
                ("proton", "V4", 4),
                ("proton", "V4", 5),
                ("proton", "V4", 6),
                ("proton", "V4", 7),
                ("proton", "V4", 8),
                ("proton", "V4", 9),
                ("proton", "V4", 10),
                #("proton", "V4", 11),
            ],
        ]
        #p["sequence"] = [[("proton", "twopt"), ("proton", "fhA3"), ("proton", "fhV4")]]
        #p["sequence"] = [[("proton", "twopt")]]
        # trange for seq3pt is for insertion-sink separation time
        p["trange"] = dict()
        t = p["trange"]
        t[("proton", "twopt")] = [3, 10]
        t[("proton", "smA3")] = [3, 10]
        t[("proton", "smV4")] = [3, 10]
        t[("proton", "fhA3")] = [2, 10]
        t[("proton", "fhV4")] = [2, 10]
        t[("proton", "A3", 3)] = [1, 2]
        t[("proton", "A3", 4)] = [2, 2]
        t[("proton", "A3", 5)] = [2, 3]
        t[("proton", "A3", 6)] = [2, 4]
        t[("proton", "A3", 7)] = [2, 5]
        t[("proton", "A3", 8)] = [2, 6]
        t[("proton", "A3", 9)] = [2, 7]
        t[("proton", "A3", 10)] = [2, 8]
        t[("proton", "A3", 11)] = [2, 9]
        t[("proton", "V4", 3)] = [1, 2]
        t[("proton", "V4", 4)] = [2, 2]
        t[("proton", "V4", 5)] = [2, 3]
        t[("proton", "V4", 6)] = [2, 4]
        t[("proton", "V4", 7)] = [2, 5]
        t[("proton", "V4", 8)] = [2, 6]
        t[("proton", "V4", 9)] = [2, 7]
        t[("proton", "V4", 10)] = [2, 8]
        t[("proton", "V4", 11)] = [2, 9]

    # priors
    # use gv.BufferDict() for log-normal priors
    # C2pt = z^2 exp(-e * t)
    import gvar as gv

    parameters["priors"] = gv.BufferDict()
    p = parameters["priors"]
    if ensemble == "a09m310":
        p["e0"] = gv.gvar(0.5, 0.05)
        p["log(e1)"] = gv.gvar(-0.86, 0.64)
        p["log(e2)"] = gv.gvar(-0.86, 0.64)
        p["log(e3)"] = gv.gvar(-0.86, 0.64)
        p["log(e4)"] = gv.gvar(-0.86, 0.64)
        p["log(e5)"] = gv.gvar(-0.86, 0.64)
        p["log(e6)"] = gv.gvar(-0.86, 0.64)
        p["log(e7)"] = gv.gvar(-0.86, 0.64)
        p["log(e8)"] = gv.gvar(-0.86, 0.64)
        p["log(e9)"] = gv.gvar(-0.86, 0.64)

        p["z0"] = gv.gvar(3.5e-4, 3.5e-5)
        p["z1"] = gv.gvar(0.0, 7e-4)
        p["z2"] = gv.gvar(0.0, 7e-4)
        p["z3"] = gv.gvar(0.0, 7e-4)
        p["z4"] = gv.gvar(0.0, 7e-4)
        p["z5"] = gv.gvar(0.0, 7e-4)
        p["z6"] = gv.gvar(0.0, 7e-4)
        p["z7"] = gv.gvar(0.0, 7e-4)
        p["z8"] = gv.gvar(0.0, 7e-4)
        p["z9"] = gv.gvar(0.0, 7e-4)

        p["a00"] = gv.gvar(1.2, 0.2)
        p["a10"] = gv.gvar(0.0, 1.0)
        p["a11"] = gv.gvar(0.0, 1.0)
        p["a20"] = gv.gvar(0.0, 1.0)
        p["a21"] = gv.gvar(0.0, 1.0)
        p["a22"] = gv.gvar(0.0, 1.0)
        p["a30"] = gv.gvar(0.0, 1.0)
        p["a31"] = gv.gvar(0.0, 1.0)
        p["a32"] = gv.gvar(0.0, 1.0)
        p["a33"] = gv.gvar(0.0, 1.0)
        p["a40"] = gv.gvar(0.0, 1.0)
        p["a41"] = gv.gvar(0.0, 1.0)
        p["a42"] = gv.gvar(0.0, 1.0)
        p["a43"] = gv.gvar(0.0, 1.0)
        p["a44"] = gv.gvar(0.0, 1.0)
        p["a50"] = gv.gvar(0.0, 1.0)
        p["a51"] = gv.gvar(0.0, 1.0)
        p["a52"] = gv.gvar(0.0, 1.0)
        p["a53"] = gv.gvar(0.0, 1.0)
        p["a54"] = gv.gvar(0.0, 1.0)
        p["a55"] = gv.gvar(0.0, 1.0)

        p["v00"] = gv.gvar(1.0, 0.2)
        p["v10"] = gv.gvar(0.0, 1.0)
        p["v11"] = gv.gvar(0.0, 1.0)
        p["v20"] = gv.gvar(0.0, 1.0)
        p["v21"] = gv.gvar(0.0, 1.0)
        p["v22"] = gv.gvar(0.0, 1.0)
        p["v30"] = gv.gvar(0.0, 1.0)
        p["v31"] = gv.gvar(0.0, 1.0)
        p["v32"] = gv.gvar(0.0, 1.0)
        p["v33"] = gv.gvar(0.0, 1.0)
        p["v40"] = gv.gvar(0.0, 1.0)
        p["v41"] = gv.gvar(0.0, 1.0)
        p["v42"] = gv.gvar(0.0, 1.0)
        p["v43"] = gv.gvar(0.0, 1.0)
        p["v44"] = gv.gvar(0.0, 1.0)
        p["v50"] = gv.gvar(0.0, 1.0)
        p["v51"] = gv.gvar(0.0, 1.0)
        p["v52"] = gv.gvar(0.0, 1.0)
        p["v53"] = gv.gvar(0.0, 1.0)
        p["v54"] = gv.gvar(0.0, 1.0)
        p["v55"] = gv.gvar(0.0, 1.0)
    elif ensemble == "a12m310":
        p["e0"] = gv.gvar(0.5, 0.05)
    elif ensemble == "a15m310":
        p["e0"] = gv.gvar(0.6, 0.2)
    elif ensemble == "a09m135":
        p["e0"] = gv.gvar(0.42, 0.042)
        p["log(e1)"] = gv.gvar(-2.1, 0.7)
        p["log(e2)"] = gv.gvar(-2.1, 0.7)
        p["log(e3)"] = gv.gvar(-2.1, 0.7)
        p["log(e4)"] = gv.gvar(-2.1, 0.7)
        p["log(e5)"] = gv.gvar(-2.1, 0.7)
        p["log(e6)"] = gv.gvar(-2.1, 0.7)
        p["log(e7)"] = gv.gvar(-2.1, 0.7)
        p["log(e8)"] = gv.gvar(-2.1, 0.7)
        p["log(e9)"] = gv.gvar(-2.1, 0.7)

        p["z0"] = gv.gvar(3e-4, 3e-5)
        p["z1"] = gv.gvar(0.0, 3e-4)
        p["z2"] = gv.gvar(0.0, 3e-4)
        p["z3"] = gv.gvar(0.0, 3e-4)
        p["z4"] = gv.gvar(0.0, 3e-4)
        p["z5"] = gv.gvar(0.0, 3e-4)
        p["z6"] = gv.gvar(0.0, 3e-4)
        p["z7"] = gv.gvar(0.0, 3e-4)
        p["z8"] = gv.gvar(0.0, 3e-4)
        p["z9"] = gv.gvar(0.0, 3e-4)

        p["a00"] = gv.gvar(1.2, 0.2)
        p["a10"] = gv.gvar(0.0, 1.0)
        p["a11"] = gv.gvar(0.0, 1.0)
        p["a20"] = gv.gvar(0.0, 1.0)
        p["a21"] = gv.gvar(0.0, 1.0)
        p["a22"] = gv.gvar(0.0, 1.0)
        p["a30"] = gv.gvar(0.0, 1.0)
        p["a31"] = gv.gvar(0.0, 1.0)
        p["a32"] = gv.gvar(0.0, 1.0)
        p["a33"] = gv.gvar(0.0, 1.0)
        p["a40"] = gv.gvar(0.0, 1.0)
        p["a41"] = gv.gvar(0.0, 1.0)
        p["a42"] = gv.gvar(0.0, 1.0)
        p["a43"] = gv.gvar(0.0, 1.0)
        p["a44"] = gv.gvar(0.0, 1.0)
        p["a50"] = gv.gvar(0.0, 1.0)
        p["a51"] = gv.gvar(0.0, 1.0)
        p["a52"] = gv.gvar(0.0, 1.0)
        p["a53"] = gv.gvar(0.0, 1.0)
        p["a54"] = gv.gvar(0.0, 1.0)
        p["a55"] = gv.gvar(0.0, 1.0)
        p["a60"] = gv.gvar(0.0, 1.0)
        p["a61"] = gv.gvar(0.0, 1.0)
        p["a62"] = gv.gvar(0.0, 1.0)
        p["a63"] = gv.gvar(0.0, 1.0)
        p["a64"] = gv.gvar(0.0, 1.0)
        p["a65"] = gv.gvar(0.0, 1.0)
        p["a66"] = gv.gvar(0.0, 1.0)
        p["a70"] = gv.gvar(0.0, 1.0)
        p["a71"] = gv.gvar(0.0, 1.0)
        p["a72"] = gv.gvar(0.0, 1.0)
        p["a73"] = gv.gvar(0.0, 1.0)
        p["a74"] = gv.gvar(0.0, 1.0)
        p["a75"] = gv.gvar(0.0, 1.0)
        p["a76"] = gv.gvar(0.0, 1.0)
        p["a77"] = gv.gvar(0.0, 1.0)

        p["v00"] = gv.gvar(1.0, 0.2)
        p["v10"] = gv.gvar(0.0, 1.0)
        p["v11"] = gv.gvar(0.0, 1.0)
        p["v20"] = gv.gvar(0.0, 1.0)
        p["v21"] = gv.gvar(0.0, 1.0)
        p["v22"] = gv.gvar(0.0, 1.0)
        p["v30"] = gv.gvar(0.0, 1.0)
        p["v31"] = gv.gvar(0.0, 1.0)
        p["v32"] = gv.gvar(0.0, 1.0)
        p["v33"] = gv.gvar(0.0, 1.0)
        p["v40"] = gv.gvar(0.0, 1.0)
        p["v41"] = gv.gvar(0.0, 1.0)
        p["v42"] = gv.gvar(0.0, 1.0)
        p["v43"] = gv.gvar(0.0, 1.0)
        p["v44"] = gv.gvar(0.0, 1.0)
        p["v50"] = gv.gvar(0.0, 1.0)
        p["v51"] = gv.gvar(0.0, 1.0)
        p["v52"] = gv.gvar(0.0, 1.0)
        p["v53"] = gv.gvar(0.0, 1.0)
        p["v54"] = gv.gvar(0.0, 1.0)
        p["v55"] = gv.gvar(0.0, 1.0)
        p["v60"] = gv.gvar(0.0, 1.0)
        p["v61"] = gv.gvar(0.0, 1.0)
        p["v62"] = gv.gvar(0.0, 1.0)
        p["v63"] = gv.gvar(0.0, 1.0)
        p["v64"] = gv.gvar(0.0, 1.0)
        p["v65"] = gv.gvar(0.0, 1.0)
        p["v66"] = gv.gvar(0.0, 1.0)
        p["v70"] = gv.gvar(0.0, 1.0)
        p["v71"] = gv.gvar(0.0, 1.0)
        p["v72"] = gv.gvar(0.0, 1.0)
        p["v73"] = gv.gvar(0.0, 1.0)
        p["v74"] = gv.gvar(0.0, 1.0)
        p["v75"] = gv.gvar(0.0, 1.0)
        p["v76"] = gv.gvar(0.0, 1.0)
        p["v77"] = gv.gvar(0.0, 1.0)
    elif ensemble == "a12m135":
        p["e0"] = gv.gvar(0.58, 0.058)
        p["log(e1)"] = gv.gvar(-1.8, 0.7)
        p["log(e2)"] = gv.gvar(-1.8, 0.7)
        p["log(e3)"] = gv.gvar(-1.8, 0.7)
        p["log(e4)"] = gv.gvar(-1.8, 0.7)
        p["log(e5)"] = gv.gvar(-1.8, 0.7)
        p["log(e6)"] = gv.gvar(-1.8, 0.7)
        p["log(e7)"] = gv.gvar(-1.8, 0.7)
        p["log(e8)"] = gv.gvar(-1.8, 0.7)
        p["log(e9)"] = gv.gvar(-1.8, 0.7)

        p["z0"] = gv.gvar(8e-4, 8e-5)
        p["z1"] = gv.gvar(4e-4, 4e-4)
        p["z2"] = gv.gvar(4e-4, 4e-4)
        p["z3"] = gv.gvar(4e-4, 4e-4)
        p["z4"] = gv.gvar(4e-4, 4e-4)
        p["z5"] = gv.gvar(4e-4, 4e-4)
        p["z6"] = gv.gvar(4e-4, 4e-4)
        p["z7"] = gv.gvar(4e-4, 4e-4)
        p["z8"] = gv.gvar(4e-4, 4e-4)
        p["z9"] = gv.gvar(4e-4, 4e-4)

        p["a00"] = gv.gvar(1.2, 0.2)
        p["a10"] = gv.gvar(0.0, 1.0)
        p["a11"] = gv.gvar(0.0, 1.0)
        p["a20"] = gv.gvar(0.0, 1.0)
        p["a21"] = gv.gvar(0.0, 1.0)
        p["a22"] = gv.gvar(0.0, 1.0)
        p["a30"] = gv.gvar(0.0, 1.0)
        p["a31"] = gv.gvar(0.0, 1.0)
        p["a32"] = gv.gvar(0.0, 1.0)
        p["a33"] = gv.gvar(0.0, 1.0)
        p["a40"] = gv.gvar(0.0, 1.0)
        p["a41"] = gv.gvar(0.0, 1.0)
        p["a42"] = gv.gvar(0.0, 1.0)
        p["a43"] = gv.gvar(0.0, 1.0)
        p["a44"] = gv.gvar(0.0, 1.0)
        p["a50"] = gv.gvar(0.0, 1.0)
        p["a51"] = gv.gvar(0.0, 1.0)
        p["a52"] = gv.gvar(0.0, 1.0)
        p["a53"] = gv.gvar(0.0, 1.0)
        p["a54"] = gv.gvar(0.0, 1.0)
        p["a55"] = gv.gvar(0.0, 1.0)
        p["a60"] = gv.gvar(0.0, 1.0)
        p["a61"] = gv.gvar(0.0, 1.0)
        p["a62"] = gv.gvar(0.0, 1.0)
        p["a63"] = gv.gvar(0.0, 1.0)
        p["a64"] = gv.gvar(0.0, 1.0)
        p["a65"] = gv.gvar(0.0, 1.0)
        p["a66"] = gv.gvar(0.0, 1.0)
        p["a70"] = gv.gvar(0.0, 1.0)
        p["a71"] = gv.gvar(0.0, 1.0)
        p["a72"] = gv.gvar(0.0, 1.0)
        p["a73"] = gv.gvar(0.0, 1.0)
        p["a74"] = gv.gvar(0.0, 1.0)
        p["a75"] = gv.gvar(0.0, 1.0)
        p["a76"] = gv.gvar(0.0, 1.0)
        p["a77"] = gv.gvar(0.0, 1.0)

        p["v00"] = gv.gvar(1.0, 0.2)
        p["v10"] = gv.gvar(0.0, 1.0)
        p["v11"] = gv.gvar(0.0, 1.0)
        p["v20"] = gv.gvar(0.0, 1.0)
        p["v21"] = gv.gvar(0.0, 1.0)
        p["v22"] = gv.gvar(0.0, 1.0)
        p["v30"] = gv.gvar(0.0, 1.0)
        p["v31"] = gv.gvar(0.0, 1.0)
        p["v32"] = gv.gvar(0.0, 1.0)
        p["v33"] = gv.gvar(0.0, 1.0)
        p["v40"] = gv.gvar(0.0, 1.0)
        p["v41"] = gv.gvar(0.0, 1.0)
        p["v42"] = gv.gvar(0.0, 1.0)
        p["v43"] = gv.gvar(0.0, 1.0)
        p["v44"] = gv.gvar(0.0, 1.0)
        p["v50"] = gv.gvar(0.0, 1.0)
        p["v51"] = gv.gvar(0.0, 1.0)
        p["v52"] = gv.gvar(0.0, 1.0)
        p["v53"] = gv.gvar(0.0, 1.0)
        p["v54"] = gv.gvar(0.0, 1.0)
        p["v55"] = gv.gvar(0.0, 1.0)
        p["v60"] = gv.gvar(0.0, 1.0)
        p["v61"] = gv.gvar(0.0, 1.0)
        p["v62"] = gv.gvar(0.0, 1.0)
        p["v63"] = gv.gvar(0.0, 1.0)
        p["v64"] = gv.gvar(0.0, 1.0)
        p["v65"] = gv.gvar(0.0, 1.0)
        p["v66"] = gv.gvar(0.0, 1.0)
        p["v70"] = gv.gvar(0.0, 1.0)
        p["v71"] = gv.gvar(0.0, 1.0)
        p["v72"] = gv.gvar(0.0, 1.0)
        p["v73"] = gv.gvar(0.0, 1.0)
        p["v74"] = gv.gvar(0.0, 1.0)
        p["v75"] = gv.gvar(0.0, 1.0)
        p["v76"] = gv.gvar(0.0, 1.0)
        p["v77"] = gv.gvar(0.0, 1.0)

        p["c0"] = gv.gvar(-1.5E-6, 1.5E-6)
        p["c1"] = gv.gvar(0.0, 1.5E-6)
        p["c2"] = gv.gvar(0.0, 1.5E-6)
        p["c3"] = gv.gvar(0.0, 1.5E-6)
        p["c4"] = gv.gvar(0.0, 1.5E-6)
        p["c5"] = gv.gvar(0.0, 1.5E-6)
        p["c6"] = gv.gvar(0.0, 1.5E-6)
        p["c7"] = gv.gvar(0.0, 1.5E-6)

        p["d0"] = gv.gvar(1.3E-6, 1.3E-6)
        p["d1"] = gv.gvar(0.0, 1.3E-6)
        p["d2"] = gv.gvar(0.0, 1.3E-6)
        p["d3"] = gv.gvar(0.0, 1.3E-6)
        p["d4"] = gv.gvar(0.0, 1.3E-6)
        p["d5"] = gv.gvar(0.0, 1.3E-6)
        p["d6"] = gv.gvar(0.0, 1.3E-6)
        p["d7"] = gv.gvar(0.0, 1.3E-6)
    return parameters
