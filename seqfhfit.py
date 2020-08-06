import numpy as np
import pandas as pd
import gvar as gv
import lsqfit
import os
import h5py
import matplotlib.pyplot as plt

import h5file_parameters as fparams


class Fitter:
    def __init__(self, fparams):
        self.strategy = fparams["strategy"]
        self.init_priors = fparams["priors"]
        self.ensemble = fparams["ensemble"]
        self.a_fm = float("0.%s" % fparams["h5filename"][1:3])
        self.mpi = float(fparams["h5filename"].split("m")[1][:3])
        self.tsep = fparams["tsep"]
        self.mval = fparams["mval"]
        self.h5filename = fparams["h5filename"]
        self.h5prepath = fparams["h5prepath"]
        self.data = self.read_data()
        self.current = np.unique([key[1] for key in self.data.keys() if len(key) == 3])
        self.nt = len(self.data[("proton", "twopt")])

    def old_fh(self):
        fhpath = f"{self.h5prepath}/fh/ml{str(self.mval).replace('.','p')}"
        fhfile = self.h5filename.replace("avg", "fh").replace("srcs", "Srcs")
        fhdata = h5py.File(f"./data/{fhfile}", "r")

        ddat = dict()
        for op in ["A3", "V4"]:
            for iso in ["UU", "DD"]:
                psavg = 0
                for parity in ["", "_np"]:
                    for spin in ["spin_up", "spin_dn"]:
                        corr = f"fh_proton{parity}_{op}_{iso}/{spin}"
                        dat = np.squeeze(fhdata[f"{fhpath}/{corr}"][:, :, 0])
                        if parity == "_np":
                            dat = np.roll(dat[:, ::-1], 1, axis=1)
                        if op == "A3":
                            dat = dat.imag
                        else:
                            dat = dat.real
                        if spin == "spin_up":
                            savg = 0.5 * dat
                        else:
                            if op == "A3":
                                savg -= 0.5 * dat
                            else:
                                savg += 0.5 * dat
                    psavg += 0.5 * savg
                if iso == "UU":
                    ipsavg = psavg
                else:
                    ipsavg -= psavg
            ddat[("proton", f"fh{op}")] = ipsavg
        return ddat

    def read_data(self):
        def data_pivot(arg):
            corr_pivot = arg.pivot(index="cfg", columns="t", values="corr").values
            return corr_pivot

        picklefile = "./data/%s.pickle" % self.ensemble
        if os.path.isfile(picklefile):
            corr_df = gv.load(picklefile)
            return corr_df
        else:
            import re
            from nucleon_elastic_ff.data.h5io import get_dsets

            # averaging three point
            pattern = "(?P<parity>proton|proton\_np)"
            pattern += "_(?P<isospin>DD|UU)"
            pattern += "_(?P<spin>dn_dn|up_up)"
            pattern += "_tsep_[\-]*(?P<tsep>[0-9]+)"
            pattern += ".*(?P<current>A3|V4).*cfgs\_srcs"

            columns = [
                "nucleon",
                "current",
                "tsep",
                "cfg",
                "t",
                "isospin",
                "parity",
                "spin",
                "corr",
            ]
            data_frames = []

            with h5py.File("./data/%s" % self.h5filename, "r") as h5f:
                dsets = get_dsets(h5f)

                for key, dset in dsets.items():
                    match = re.search(pattern, key)
                    if match:
                        info = match.groupdict()

                        nucleon_parity = info.pop("parity").split("_")
                        info["nucleon"] = nucleon_parity[0]

                        spin = info.pop("spin")
                        info["spin"] = (
                            -1 if (spin == "dn_dn" and info["current"] in ["A3"]) else 1
                        )

                        info["parity"] = -1 if len(nucleon_parity) == 2 else 1

                        isospin = info.pop("isospin")
                        info["isospin"] = 1 if isospin == "UU" else -1

                        current_key = key.replace("cfgs_srcs", "local_curr")
                        curr_dset = h5f[current_key]

                        cfgs = dset[:, 0]
                        corr = (
                            curr_dset[()].real
                            if info["current"] in ["V4"]
                            else curr_dset[()].imag
                        )
                        ts = range(corr.shape[-1])

                        tmp_df = (
                            pd.DataFrame(index=cfgs, columns=ts, data=corr)
                            .unstack()
                            .reset_index()
                            .rename(
                                columns={"level_0": "t", "level_1": "cfg", 0: "corr"}
                            )
                        )
                        for key, val in info.items():
                            tmp_df[key] = val
                        data_frames.append(tmp_df.astype({"tsep": int}))

            df = (
                pd.concat(data_frames, ignore_index=True)
                .reindex(columns, axis=1)
                .sort_values(columns)
                .reset_index(drop=True)
            )

            # spin average
            tmp = df.copy()
            tmp["corr"] *= tmp["spin"]
            spin_avg_df = tmp.groupby(
                ["nucleon", "current", "tsep", "cfg", "t", "isospin", "parity"],
                as_index=False,
            )["corr"].mean()

            # parity average
            tmp = spin_avg_df.copy()
            tmp["corr"] *= tmp["parity"]
            spin_parity_avg_df = tmp.groupby(
                ["nucleon", "current", "tsep", "cfg", "t", "isospin"], as_index=False
            )[["corr"]].mean()

            # isovector combination
            tmp = spin_parity_avg_df.copy()
            tmp["corr"] *= tmp["isospin"]
            isospin_spin_parity_avg_df = tmp.groupby(
                ["nucleon", "current", "tsep", "cfg", "t"], as_index=False
            )["corr"].sum()

            # pivot to dictionary with key = (nucleon, current, tsep), data = [cfg, t]
            group = isospin_spin_parity_avg_df.groupby(["nucleon", "current", "tsep"])
            corr_df = group.apply(data_pivot).to_dict()

            # averaging two point
            pattern = "(?P<parity>proton|proton\_np)"
            pattern += "/.*/(?P<spin>spin\_dn|spin\_up)"

            columns = ["nucleon", "cfg", "t", "parity", "spin", "corr"]
            data_frames = []

            with h5py.File("./data/%s" % self.h5filename, "r") as h5f:
                dsets = get_dsets(h5f)

                for key, dset in dsets.items():
                    match = re.search(pattern, key)
                    if match:
                        info = match.groupdict()

                        nucleon_parity = info.pop("parity").split("_")
                        info["nucleon"] = nucleon_parity[0]
                        info["parity"] = -1 if len(nucleon_parity) == 2 else 1

                        cfg_key = key.replace("spin_dn", "cfgs_srcs").replace(
                            "spin_up", "cfgs_srcs"
                        )
                        cfg_dset = h5f[cfg_key]

                        cfgs = cfg_dset[:, 0]
                        corr = dset[()][:, :, 0, 0].real
                        ts = range(corr.shape[-1])

                        tmp_df = (
                            pd.DataFrame(index=cfgs, columns=ts, data=corr)
                            .unstack()
                            .reset_index()
                            .rename(
                                columns={"level_0": "t", "level_1": "cfg", 0: "corr"}
                            )
                        )
                        for key, val in info.items():
                            tmp_df[key] = val
                        data_frames.append(tmp_df)
            df = (
                pd.concat(data_frames, ignore_index=True)
                .reindex(columns, axis=1)
                .sort_values(columns)
                .reset_index(drop=True)
            )
            # spin average
            spin_avg_df = df.groupby(["nucleon", "cfg", "t", "parity"], as_index=False)[
                "corr"
            ].mean()

            # parity average (negative sign already accounted for in run scripts)
            spin_parity_avg_df = spin_avg_df.groupby(
                ["nucleon", "cfg", "t"], as_index=False
            )[["corr"]].mean()
            spin_parity_avg_df["current"] = "twopt"

            # pivot to dictionary with key = (nucleon), data = [cfg, t]
            group = spin_parity_avg_df.groupby(["nucleon", "current"])
            corr_df.update(group.apply(data_pivot).to_dict())

            # read old style FH correlator
            if True:
                ddat = self.old_fh()
                corr_df.update(ddat)

            # make gvar
            corr_dict = gv.dataset.avg_data(corr_df)

            if True:
                # make old style fh correlator
                for op in ["A3", "V4"]:
                    R = (
                        corr_dict[("proton", f"fh{op}")]
                        / corr_dict[("proton", "twopt")]
                    )
                    fh = np.roll(R, -1) - R
                    corr_dict.update({("proton", f"fh{op}"): fh})

            # make sum correlator
            current = np.unique([key[1] for key in corr_dict.keys() if len(key) == 3])
            summation = {
                cur: np.pad(
                    [sum(corr_dict["proton", cur, t][1 : t - 1]) for t in self.tsep],
                    (
                        self.tsep[0],
                        len(corr_dict[("proton", "twopt")]) - self.tsep[-1] - 1,
                    ),
                    "constant",
                    constant_values=(0, 0),
                )
                / corr_dict[("proton", "twopt")]
                for cur in current
            }
            corr_dict.update(
                {
                    ("proton", "sm%s" % cur): np.roll(summation[cur], -1)
                    - summation[cur]
                    for cur in current
                }
            )

            # write to pickle
            gv.dump(corr_dict, picklefile)
            return corr_dict

    def plot_data(self):
        # plot two point
        twopt = self.data[("proton", "twopt")]
        effective_mass = np.log(np.roll(twopt, 1) / twopt)
        fig = plt.figure("effective mass", figsize=(7, 4))
        ax = plt.axes([0.15, 0.15, 0.8, 0.8])
        ax.errorbar(
            x=np.arange(len(effective_mass)),
            y=[i.mean for i in effective_mass],
            yerr=[i.sdev for i in effective_mass],
            marker="o",
            mfc="none",
            ls="none",
        )
        plt.draw()
        plt.show()

        effz = np.sqrt(twopt * np.exp(effective_mass * np.arange(len(effective_mass))))
        effective_z = effz[: len(effective_mass) // 4]
        fig = plt.figure("effective z", figsize=(7, 4))
        ax = plt.axes([0.15, 0.15, 0.8, 0.8])
        ax.errorbar(
            x=np.arange(len(effective_z)),
            y=[i.mean for i in effective_z],
            yerr=[i.sdev for i in effective_z],
            marker="o",
            mfc="none",
            ls="none",
        )
        plt.draw()
        plt.show()

        # plot three point
        for cur in self.current:
            fig = plt.figure("three point %s" % cur, figsize=(7, 4))
            ax = plt.axes([0.15, 0.15, 0.8, 0.8])
            for t in self.tsep:
                threept = (
                    self.data[("proton", cur, t)] / self.data[("proton", "twopt")][t]
                )
                mask = np.arange(t + 1)
                ax.errorbar(
                    x=mask - t * 0.5,
                    y=np.array([i.mean for i in threept])[mask],
                    yerr=np.array([i.sdev for i in threept])[mask],
                    marker="o",
                    mfc="none",
                    ls="none",
                )
            plt.draw()
            plt.show()

        # plot sum and fh
        for cur in self.current:
            fh = self.data[("proton", "fh%s" % cur)]
            # sum = self.data[("proton", "sm%s" % cur)]
            fig = plt.figure("sm %s" % cur, figsize=(7, 4))
            ax = plt.axes([0.15, 0.15, 0.8, 0.8])
            # ax.errorbar(
            #    x=np.arange(len(sum)),
            #    y=[i.mean for i in sum],
            #    yerr=[i.sdev for i in sum],
            #    marker="o",
            #    mfc="none",
            #    ls="none",
            #    label="sum",
            # )
            ax.errorbar(
                x=np.arange(len(fh)),
                y=[i.mean for i in fh],
                yerr=[i.sdev for i in fh],
                marker="o",
                mfc="none",
                ls="none",
                label="fh",
            )
            plt.legend()
            plt.draw()
            plt.show()

        oldfh = gv.dataset.avg_data(self.old_fh())
        for corr in oldfh:
            # print contact term strength
            E0p = self.init_priors["e0"]
            print(f"{corr[1]} contact strength: {oldfh[corr][1]*np.exp(E0p)}")

    def priors(self):
        def uncorrelate_prior(p):
            return gv.gvar(p.mean, p.sdev)
            # return p

        prior_key = dict()
        prior_key["twopt"] = ["z", "e"]
        prior_key["A3"] = ["z", "e", "a"]
        prior_key["V4"] = ["z", "e", "v"]
        prior_key["fhA3"] = ["z", "e", "a"]
        prior_key["fhV4"] = ["z", "e", "v"]
        prior_key["smA3"] = ["z", "e", "a"]
        prior_key["smV4"] = ["z", "e", "v"]
        priors = gv.BufferDict()
        flat_sequence = []
        for corr_set in self.strategy["sequence"]:
            for corr in corr_set:
                flat_sequence.append(corr[1])
        for corr in flat_sequence:
            if corr in ["fhA4", "fhV4", "smA3", "smV4"]:
                nstate = self.strategy["fstates"]
            else:
                nstate = self.strategy["nstates"]
            for p in self.init_priors:
                if p.replace("log(", "")[0] in prior_key[corr] and int(
                    p.replace("log(", "")[1]
                ) < int(nstate):
                    priors[p] = self.init_priors[p]
        # add FH trashcan if applicable
        for tag in ["sm", "fh"]:
            if (
                set((f"{tag}A3", f"{tag}V4")).intersection(set(flat_sequence))
                and self.strategy["fstates"] < self.strategy["nstates"]
            ):
                for cur in self.current:
                    priors[
                        f"log({tag}e%s)" % (self.strategy["fstates"] - 1)
                    ] = uncorrelate_prior(
                        self.init_priors["log(e%s)" % (self.strategy["fstates"] - 1)]
                    )
                    priors[
                        f"{tag}z%s" % (self.strategy["fstates"] - 1)
                    ] = uncorrelate_prior(
                        self.init_priors["z%s" % (self.strategy["fstates"] - 1)]
                    )
                    for ns in range(self.strategy["fstates"]):
                        priors[
                            f"{tag}%s%s%s"
                            % (cur[0].lower(), self.strategy["fstates"] - 1, ns)
                        ] = uncorrelate_prior(
                            self.init_priors[
                                "%s%s%s"
                                % (cur[0].lower(), self.strategy["fstates"] - 1, ns)
                            ]
                        )
        # FH contact and outside terms
        for tag in ["fhA3", "fhV4"]:
            if tag in flat_sequence:
                if tag == "fhA3":
                    label = "c"
                else:
                    label = "d"
                for ns in range(self.strategy["fstates"]):
                    priors[f"{label}{ns}"] = self.init_priors[f"{label}{ns}"]
        return priors

    def prepare_data(self, subset):
        trange = self.strategy["trange"]
        X = dict()
        y = dict()
        for corr in subset:
            X[corr] = np.arange(trange[corr][0], trange[corr][1] + 1)
            y[corr] = self.data[corr][X[corr]]
        return X, y

    def fit(self):
        p = self.priors()
        init = dict()
        init["nt"] = self.nt
        init["strategy"] = self.strategy
        self.model = FitFunction(init)
        for subset in self.strategy["sequence"]:
            X, y = self.prepare_data(subset)
            result = lsqfit.nonlinear_fit(
                data=(X, y), prior=p, fcn=self.model, maxit=10000000
            )
            p = result.p
            print(result)
        self.result = result
        print(result.format(maxline=True))

    def plot_result(self):
        step = 10  # density for fit result
        flatcorr = []
        for si in self.strategy["sequence"]:
            for s in si:
                flatcorr.append(s)
        # twopt, A3, V4, smA3, smV4, fhA3, fhV4
        def twopt():
            corrkey = ("proton", "twopt")
            twopt = self.data[corrkey]
            dmeff = np.log(np.roll(twopt, 1) / twopt)
            x = np.linspace(0, len(twopt), len(twopt) * step + 1)
            y = list(self.model({corrkey: x}, self.result.p).values())[0]
            fmeff = np.log(np.roll(y, step) / y)
            ym = np.array([i.mean for i in fmeff])
            ys = np.array([i.sdev for i in fmeff])

            plt.figure("effective mass", figsize=(7, 4))
            ax = plt.axes([0.15, 0.15, 0.8, 0.8])
            ax.errorbar(
                x=range(len(twopt)),
                y=[i.mean for i in dmeff],
                yerr=[i.sdev for i in dmeff],
                ls="none",
            )
            ax.fill_between(x=x, y1=ym - ys, y2=ym + ys, alpha=0.5)
            ax.set_ylim([0.5, 0.7])
            ax.set_xlim([2, 13])
            plt.draw()
            plt.show()

        def fh(corrkey):
            self.model.corr = corrkey
            fh = self.data[corrkey]
            x = np.linspace(0, len(fh), len(fh) * step + 1)
            numerator = self.model.summation(x, self.result.p, "fh")
            denominator = self.model.twopoint(x, self.result.p)
            ratio = numerator / denominator
            y = np.roll(ratio, -step) - ratio
            ym = np.array([i.mean for i in y])
            ys = np.array([i.sdev for i in y])

            plt.figure(corrkey[1], figsize=(7, 4))
            ax = plt.axes([0.15, 0.15, 0.8, 0.8])
            ax.errorbar(
                x=range(len(fh)),
                y=[i.mean for i in fh],
                yerr=[i.sdev for i in fh],
                ls="none",
            )
            ax.fill_between(x=x, y1=ym - ys, y2=ym + ys, alpha=0.5)
            if corrkey[1] == "fhA3":
                ax.set_ylim([1.15, 1.4])
            elif corrkey[1] == "fhV4":
                ax.set_ylim([0.98, 1.1])
            ax.set_xlim([2, 13])
            plt.draw()
            plt.show()

        def seq(glist):
            xi = {key: np.array(range(0, key[2] + 1)) for key in glist}
            di = {
                key: self.data[key] / self.data[("proton", "twopt")][key[2]]
                for key in glist
            }
            fx = {key: np.linspace(0, key[2], key[2] * step + 1) for key in glist}
            twopt = list(
                self.model(
                    {("proton", "twopt"): np.array(range(self.nt))}, self.result.p
                ).values()
            )[0]
            f3pt = self.model(fx, self.result.p)
            y = {key: f3pt[key] / twopt[key[2]] for key in f3pt}

            plt.figure(glist[0][1], figsize=(7, 4))
            ax = plt.axes([0.15, 0.15, 0.8, 0.8])
            for key in xi:
                ax.errorbar(
                    x=xi[key] - key[2] * 0.5,
                    y=[i.mean for i in di[key]],
                    yerr=[i.sdev for i in di[key]],
                    ls="none",
                )
            for key in fx:
                yr = y[key]
                ym = np.array([i.mean for i in yr])
                ys = np.array([i.sdev for i in yr])
                ax.fill_between(
                    x=fx[key] - key[2] * 0.5, y1=ym - ys, y2=ym + ys, alpha=0.5
                )
            plt.draw()
            plt.show()

        a3list = []
        v4list = []
        for corrkey in flatcorr:
            corr = corrkey[1]
            if corr == "twopt":
                twopt()
            elif corr == "A3":
                a3list.append(corrkey)
            elif corr == "V4":
                v4list.append(corrkey)
            elif corr == "fhA3":
                fh(("proton", "fhA3"))
            elif corr == "fhV4":
                fh(("proton", "fhV4"))
        if a3list:
            seq(a3list)
        if v4list:
            seq(v4list)


class FitFunction:
    def __init__(self, init):
        self.nt = init["nt"]
        self.strategy = init["strategy"]

    def __call__(self, X, p):
        r = dict()
        for corr in X:
            self.corr = corr
            if corr[1] in ["twopt"]:
                r[corr] = self.twopoint(X[corr], p)
            elif corr[1] in ["smA3", "smV4"]:
                r[corr] = self.fh(X[corr], p, "sm")
            elif corr[1] in ["A3", "V4"]:
                r[corr] = self.threepoint(X[corr], p)
            elif corr[1] in ["fhA3", "fhV4"]:
                r[corr] = self.fh(X[corr], p, "fh")
        return r

    def twopoint(self, Xi, p):
        r = 0
        for ni in range(self.strategy["nstates"]):
            En = sum([p["e%s" % ntower] for ntower in range(ni + 1)])
            r += p["z%s" % ni] ** 2 * np.exp(-En * Xi)
        return r

    def summation(self, Xi, p, tag):
        r = 0
        for ni in range(self.strategy["fstates"]):
            for mi in range(self.strategy["fstates"]):
                sort_idx = np.sort([mi, ni])
                zn = p["z%s" % ni]
                zm = p["z%s" % mi]
                En_list = [p["e%s" % ntower] for ntower in range(ni + 1)]
                Em_list = [p["e%s" % mtower] for mtower in range(mi + 1)]
                Omn = p["%s%s%s" % (self.corr[1][2].lower(), sort_idx[1], sort_idx[0])]
                if (
                    self.strategy["fstates"] < self.strategy["nstates"]
                    and ni == self.strategy["fstates"] - 1
                ):
                    zn = p[f"{tag}z%s" % ni]
                    En_list[-1] = p[f"{tag}e%s" % ni]
                    Omn = p[
                        f"{tag}%s%s%s"
                        % (self.corr[1][2].lower(), sort_idx[1], sort_idx[0])
                    ]
                if (
                    self.strategy["fstates"] < self.strategy["nstates"]
                    and mi == self.strategy["fstates"] - 1
                ):
                    zm = p[f"{tag}z%s" % mi]
                    Em_list[-1] = p[f"{tag}e%s" % mi]
                    Omn = p[
                        f"{tag}%s%s%s"
                        % (self.corr[1][2].lower(), sort_idx[1], sort_idx[0])
                    ]
                En = sum(En_list)
                Em = sum(Em_list)
                if ni == mi:
                    r += (Xi - 1.0) * zm * Omn * zn * np.exp(-En * Xi)
                else:
                    Dmn = Em - En
                    Dnm = En - Em
                    r += (
                        zm
                        * Omn
                        * zn
                        * (
                            np.exp(-En * Xi) * np.exp(Dnm / 2.0)
                            - np.exp(-Em * Xi) * np.exp(Dmn / 2.0)
                        )
                        / (np.exp(Dmn / 2.0) - np.exp(Dnm / 2.0))
                    )
            if tag == "fh":
                if self.corr[1][2].lower() == "a":
                    label = "c"
                else:
                    label = "d"
                Dn = p[f"{label}{ni}"]
                r += Dn * np.exp(-En * Xi)
        return r

    def fh(self, Xi, p, tag):
        numerator = self.summation(np.arange(self.nt), p, tag)
        denominator = self.twopoint(np.arange(self.nt), p)
        ratio = numerator / denominator
        r = (np.roll(ratio, -1) - ratio)[Xi]
        #r = (ratio - np.roll(ratio, 1))[Xi]
        return r

    def threepoint(self, Xi, p):
        r = 0
        for ni in range(self.strategy["nstates"]):
            for mi in range(self.strategy["nstates"]):
                sort_idx = np.sort([mi, ni])
                zm = p["z%s" % mi]
                zn = p["z%s" % ni]
                En = sum([p["e%s" % ntower] for ntower in range(ni + 1)])
                Em = sum([p["e%s" % mtower] for mtower in range(mi + 1)])
                Omn = p["%s%s%s" % (self.corr[1][0].lower(), sort_idx[1], sort_idx[0])]
                r += zm * Omn * zn * np.exp(-Em * float(self.corr[2]) - (En - Em) * Xi)
        return r


if __name__ == "__main__":
    print(gv.__version__)
    ensemble = "a12m135"
    # instantiate data
    fitter = Fitter(fparams.read_parameters(ensemble))
    # print(fitter.data[("proton", "fhA3")])
    # plot data
    # fitter.plot_data()
    # fit data
    fitter.fit()
    fitter.plot_result()
    # print(fitter.to_lattedb)
