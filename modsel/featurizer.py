import pandas as pd
from warnings import filterwarnings

filterwarnings("ignore")


def gen_feats(df):

    df["nbo*"] = df["nbo1"] * df["nbo2"]
    df["nbo-"] = abs(df["nbo2"] - df["nbo1"])
    df["nbo+"] = df["nbo2"] + df["nbo1"]
    df["nbom"] = (df["nbo2"] + df["nbo1"]) / 2
    df["nbo1/2"] = df["nbo1"] / df["nbo2"]

    df["cnbo*"] = df["cnbo1"] * df["cnbo2"]
    df["cnbo-"] = abs(df["cnbo2"] - df["cnbo1"])
    df["cnbo+"] = df["cnbo2"] + df["cnbo1"]
    df["cnbom"] = (df["cnbo2"] + df["cnbo1"]) / 2
    df["cnbo1/2"] = df["cnbo1"] / df["cnbo2"]

    df["onbo*"] = df["onbo1"] * df["onbo2"]
    df["onbo-"] = abs(df["onbo2"] - df["onbo1"])
    df["onbo+"] = df["onbo2"] + df["onbo1"]
    df["onbom"] = (df["onbo2"] + df["onbo1"]) / 2
    df["onbo1/2"] = df["onbo1"] / df["onbo2"]

    df["gap"] = df["lumo"] - df["homo"]
    # df["gap*"] = df["lumo"] * df["homo"]

    df["festrada"] = df["estrada"] / df["na_an"]
    df["fwiener"] = df["wiener"] / df["na_an"]
    df["fglobal_eff"] = df["global_eff"] / df["na_an"]
    df["fbalaban"] = df["balaban"] / df["na_an"]
    df["fhosoya"] = df["hosoya"] / df["na_an"]
    df["fzagreb1"] = df["zagreb1"] / df["na_an"]
    df["fzagreb2"] = df["zagreb2"] / df["na_an"]

    df["fkier_a"] = df["kier_a"] / df["na_an"]
    df["fkier_b"] = df["kier_b"] / df["na_an"]
    df["fkier_al"] = df["kier_al"] / df["na_an"]

    df["f0k"] = df["0k"] / df["na_an"]
    df["f1k"] = df["1k"] / df["na_an"]
    df["f2k"] = df["2k"] / df["na_an"]
    df["f3k"] = df["3k"] / df["na_an"]
    df["f1ka"] = df["1ka"] / df["na_an"]
    df["f2ka"] = df["2ka"] / df["na_an"]
    df["f3ka"] = df["3ka"] / df["na_an"]
    df["f1kb"] = df["1kb"] / df["na_an"]
    df["f2kb"] = df["2kb"] / df["na_an"]
    df["f3kb"] = df["3kb"] / df["na_an"]
    df["f1kal"] = df["1kal"] / df["na_an"]
    df["f2kal"] = df["2kal"] / df["na_an"]
    df["f3kal"] = df["3kal"] / df["na_an"]

    df["f0chi"] = df["0chi"] / df["na_an"]
    df["f1chi"] = df["1chi"] / df["na_an"]
    df["f2chi"] = df["2chi"] / df["na_an"]
    df["f3chi"] = df["3chi"] / df["na_an"]
    df["f4chi"] = df["4chi"] / df["na_an"]
    df["f5chi"] = df["5chi"] / df["na_an"]
    df["f0chiv"] = df["0chiv"] / df["na_an"]
    df["f1chiv"] = df["1chiv"] / df["na_an"]
    df["f2chiv"] = df["2chiv"] / df["na_an"]
    df["f3chiv"] = df["3chiv"] / df["na_an"]
    df["f4chiv"] = df["4chiv"] / df["na_an"]
    df["f5chiv"] = df["5chiv"] / df["na_an"]

    df["d0chi"] = df["0chi"] - df["0chiv"]
    df["d1chi"] = df["1chi"] - df["1chiv"]
    df["d2chi"] = df["2chi"] - df["2chiv"]
    df["d3chi"] = df["3chi"] - df["3chiv"]
    df["d4chi"] = df["4chi"] - df["4chiv"]
    df["d5chi"] = df["5chi"] - df["5chiv"]

    df["s0chi"] = df["0chi"] + df["0chiv"]
    df["s1chi"] = df["1chi"] + df["1chiv"]
    df["s2chi"] = df["2chi"] + df["2chiv"]
    df["s3chi"] = df["3chi"] + df["3chiv"]
    df["s4chi"] = df["4chi"] + df["4chiv"]
    df["s5chi"] = df["5chi"] + df["5chiv"]

    df["fT"] = df["T"] / df["na_an"]

    df["k_ups"] = (df["1k"] * df["2k"] * df["3k"]) / df["na_an"]
    df["k_upsa"] = (df["1ka"] * df["2ka"] * df["3ka"]) / df["na_an"]
    df["k_upsb"] = (df["1kb"] * df["2kb"] * df["3kb"]) / df["na_an"]
    df["k_upsal"] = (df["1kal"] * df["2kal"] * df["3kal"]) / df["na_an"]

    df = df.copy()

    # b noh oct
    df["bonh_t"] = (
        df["bonh1"]
        + df["bonh2"]
        + df["bonh3"]
        + df["bonh4"]
        + df["bonh5"]
        + df["bonh6"]
        + df["bonh7"]
        + df["bonh8"]
    )
    df["bonhq1"] = df["bonh1"] + df["bonh8"]
    df["bonhq2"] = df["bonh2"] + df["bonh7"]
    df["bonhq3"] = df["bonh3"] + df["bonh6"]
    df["bonhq4"] = df["bonh4"] + df["bonh5"]

    df["fbonhq1"] = df["bonhq1"] / df["bonh_t"]
    df["fbonhq2"] = df["bonhq2"] / df["bonh_t"]
    df["fbonhq3"] = df["bonhq3"] / df["bonh_t"]
    df["fbonhq4"] = df["bonhq4"] / df["bonh_t"]

    # noh oct
    df["onh_t"] = (
        df["onh1"]
        + df["onh2"]
        + df["onh3"]
        + df["onh4"]
        + df["onh5"]
        + df["onh6"]
        + df["onh7"]
        + df["onh8"]
    )
    df["onhq1"] = df["onh1"] + df["onh8"]
    df["onhq2"] = df["onh2"] + df["onh7"]
    df["onhq3"] = df["onh3"] + df["onh6"]
    df["onhq4"] = df["onh4"] + df["onh5"]

    df["fonhq1"] = df["onhq1"] / df["onh_t"]
    df["fonhq2"] = df["onhq2"] / df["onh_t"]
    df["fonhq3"] = df["onhq3"] / df["onh_t"]
    df["fonhq4"] = df["onhq4"] / df["onh_t"]

    # oct
    df["o_t"] = (
        df["o1"]
        + df["o2"]
        + df["o3"]
        + df["o4"]
        + df["o5"]
        + df["o6"]
        + df["o7"]
        + df["o8"]
    )
    df["oq1"] = df["o1"] + df["o8"]
    df["oq2"] = df["o2"] + df["o7"]
    df["oq3"] = df["o3"] + df["o6"]
    df["oq4"] = df["o4"] + df["o5"]

    df["foq1"] = df["oq1"] / df["o_t"]
    df["foq2"] = df["oq2"] / df["o_t"]
    df["foq3"] = df["oq3"] / df["o_t"]
    df["foq4"] = df["oq4"] / df["o_t"]

    # Halves
    df["bonhh1"] = df["bonh1"] + df["bonh4"] + df["bonh5"] + df["bonh8"]
    df["bonhh2"] = df["bonh1"] + df["bonh2"] + df["bonh7"] + df["bonh8"]
    df["bonhh3"] = df["bonh1"] + df["bonh2"] + df["bonh3"] + df["bonh4"]
    df["bonhh4"] = df["bonh2"] + df["bonh3"] + df["bonh6"] + df["bonh7"]
    df["bonhh5"] = df["bonh3"] + df["bonh4"] + df["bonh5"] + df["bonh6"]
    df["bonhh6"] = df["bonh5"] + df["bonh6"] + df["bonh7"] + df["bonh8"]

    df["fbonhh1"] = df["bonhh1"] / df["bonh_t"]
    df["fbonhh2"] = df["bonhh2"] / df["bonh_t"]
    df["fbonhh3"] = df["bonhh3"] / df["bonh_t"]
    df["fbonhh4"] = df["bonhh4"] / df["bonh_t"]
    df["fbonhh5"] = df["bonhh5"] / df["bonh_t"]
    df["fbonhh6"] = df["bonhh6"] / df["bonh_t"]

    # Halves
    df["onhh1"] = df["onh1"] + df["onh4"] + df["onh5"] + df["onh8"]
    df["onhh2"] = df["onh1"] + df["onh2"] + df["onh7"] + df["onh8"]
    df["onhh3"] = df["onh1"] + df["onh2"] + df["onh3"] + df["onh4"]
    df["onhh4"] = df["onh2"] + df["onh3"] + df["onh6"] + df["onh7"]
    df["onhh5"] = df["onh3"] + df["onh4"] + df["onh5"] + df["onh6"]
    df["onhh6"] = df["onh5"] + df["onh6"] + df["onh7"] + df["onh8"]

    df["fonhh1"] = df["onhh1"] / df["onh_t"]
    df["fonhh2"] = df["onhh2"] / df["onh_t"]
    df["fonhh3"] = df["onhh3"] / df["onh_t"]
    df["fonhh4"] = df["onhh4"] / df["onh_t"]
    df["fonhh5"] = df["onhh5"] / df["onh_t"]
    df["fonhh6"] = df["onhh6"] / df["onh_t"]

    # Halves
    df["oh1"] = df["o1"] + df["o4"] + df["o5"] + df["o8"]
    df["oh2"] = df["o1"] + df["o2"] + df["o7"] + df["o8"]
    df["oh3"] = df["o1"] + df["o2"] + df["o3"] + df["o4"]
    df["oh4"] = df["o2"] + df["o3"] + df["o6"] + df["o7"]
    df["oh5"] = df["o3"] + df["o4"] + df["o5"] + df["o6"]
    df["oh6"] = df["o5"] + df["o6"] + df["o7"] + df["o8"]

    df["foh1"] = df["oh1"] / df["o_t"]
    df["foh2"] = df["oh2"] / df["o_t"]
    df["foh3"] = df["oh3"] / df["o_t"]
    df["foh4"] = df["oh4"] / df["o_t"]
    df["foh5"] = df["oh5"] / df["o_t"]
    df["foh6"] = df["oh6"] / df["o_t"]

    # df["ohx"] = df["oh1"] / df["oh4"]
    # df["ohy"] = df["oh2"] / df["oh5"]
    # df["ohz"] = df["oh3"] / df["oh6"]

    # df["ohx*"] = df["oh1"] * df["oh4"]
    # df["ohy*"] = df["oh2"] * df["oh5"]
    # df["ohz*"] = df["oh3"] * df["oh6"]

    # df["ons"] = df["oh6"] - df["oh3"]
    # df["onhns"] = df["onhh6"] - df["onhh3"]
    # df["bonhns"] = df["bonhh6"] - df["bonhh3"]

    df["bonhr"] = (df["bonhq1"] + df["bonhq3"]) / (df["bonhq2"] + df["bonhq4"])
    df["onhr"] = (df["onhq1"] + df["onhq3"]) / (df["onhq2"] + df["onhq4"])
    df["or"] = (df["oq1"] + df["oq3"]) / (df["oq2"] + df["oq4"])

    df["pb_vol"] = df["b_vol"] / df["t_vol"]
    df["pb_vol_nh"] = df["b_vol_nh"] / df["t_vol_nh"]

    return df.copy()


def bl_featurize(lig_feats, lig_resps=None, uid="c_smiles", task=None, target="ddg", full=False):

    df = pd.merge(lig_resps, lig_feats, on=uid)
    df_bu = df.copy()  # df back up

    for row in df_bu.iterrows():

        # if (row[1]["o1"] + row[1]["o4"] + row[1]["o5"] + row[1]["o8"]) < (
        #     row[1]["o2"] + row[1]["o3"] + row[1]["o6"] + row[1]["o7"]
        # ):
        if False:
            df.at[row[0], "bonh1"] = row[1]["bonh3"]
            df.at[row[0], "bonh2"] = row[1]["bonh4"]
            df.at[row[0], "bonh3"] = row[1]["bonh1"]
            df.at[row[0], "bonh4"] = row[1]["bonh2"]
            df.at[row[0], "bonh5"] = row[1]["bonh7"]
            df.at[row[0], "bonh6"] = row[1]["bonh8"]
            df.at[row[0], "bonh7"] = row[1]["bonh5"]
            df.at[row[0], "bonh8"] = row[1]["bonh6"]

            df.at[row[0], "onh1"] = row[1]["onh3"]
            df.at[row[0], "onh2"] = row[1]["onh4"]
            df.at[row[0], "onh3"] = row[1]["onh1"]
            df.at[row[0], "onh4"] = row[1]["onh2"]
            df.at[row[0], "onh5"] = row[1]["onh7"]
            df.at[row[0], "onh6"] = row[1]["onh8"]
            df.at[row[0], "onh7"] = row[1]["onh5"]
            df.at[row[0], "onh8"] = row[1]["onh6"]

            df.at[row[0], "o1"] = row[1]["o3"]
            df.at[row[0], "o2"] = row[1]["o4"]
            df.at[row[0], "o3"] = row[1]["o1"]
            df.at[row[0], "o4"] = row[1]["o2"]
            df.at[row[0], "o5"] = row[1]["o7"]
            df.at[row[0], "o6"] = row[1]["o8"]
            df.at[row[0], "o7"] = row[1]["o5"]
            df.at[row[0], "o8"] = row[1]["o6"]

            df.at[row[0], "nbo1"] = row[1]["nbo2"]
            df.at[row[0], "nbo2"] = row[1]["nbo1"]
            df.at[row[0], "cnbo1"] = row[1]["cnbo2"]
            df.at[row[0], "cnbo2"] = row[1]["cnbo1"]
            df.at[row[0], "onbo1"] = row[1]["onbo2"]
            df.at[row[0], "onbo2"] = row[1]["onbo1"]

        if row[1][target] < 0:
            df.at[row[0], target] = row[1][target] * -1

            df.at[row[0], "bonh1"] = row[1]["bonh4"]
            df.at[row[0], "bonh2"] = row[1]["bonh3"]
            df.at[row[0], "bonh3"] = row[1]["bonh2"]
            df.at[row[0], "bonh4"] = row[1]["bonh1"]
            df.at[row[0], "bonh5"] = row[1]["bonh8"]
            df.at[row[0], "bonh6"] = row[1]["bonh7"]
            df.at[row[0], "bonh7"] = row[1]["bonh6"]
            df.at[row[0], "bonh8"] = row[1]["bonh5"]

            df.at[row[0], "onh1"] = row[1]["onh4"]
            df.at[row[0], "onh2"] = row[1]["onh3"]
            df.at[row[0], "onh3"] = row[1]["onh2"]
            df.at[row[0], "onh4"] = row[1]["onh1"]
            df.at[row[0], "onh5"] = row[1]["onh8"]
            df.at[row[0], "onh6"] = row[1]["onh7"]
            df.at[row[0], "onh7"] = row[1]["onh6"]
            df.at[row[0], "onh8"] = row[1]["onh5"]

            df.at[row[0], "o1"] = row[1]["o4"]
            df.at[row[0], "o2"] = row[1]["o3"]
            df.at[row[0], "o3"] = row[1]["o2"]
            df.at[row[0], "o4"] = row[1]["o1"]
            df.at[row[0], "o5"] = row[1]["o8"]
            df.at[row[0], "o6"] = row[1]["o7"]
            df.at[row[0], "o7"] = row[1]["o6"]
            df.at[row[0], "o8"] = row[1]["o5"]

    df = df.query(f"set in @task")
    df = df.drop_duplicates(subset=["c_smiles"], keep="last")

    df_defrag = gen_feats(df)

    feats = [
        x
        for x in df_defrag.columns.to_list()
        if x
        not in [
            "source",
            "smiles",
            "c_smiles",
            "set",
            "name",
            "ee",
            "ddg",
            "temp",
            "class",
            "order"
        ]
    ]

    x = df_defrag[feats]
    y = df_defrag[target]

    if full:
        return df_defrag, feats
    else:
        return x, y


def bl_pool_featurize(df, symm=False, full=False):

    df = df.copy()

    if symm:
        df_bu = df.copy()  # df back up

        for row in df_bu.iterrows():
            if (row[1]["bonh1"] + row[1]["bonh8"] + row[1]["bonh3"] + row[1]["bonh6"]) / (row[1]["bonh2"] + row[1]["bonh7"] + row[1]["bonh4"] + row[1]["bonh5"]) > 1:
                df.at[row[0], "bonh1"] = row[1]["bonh4"]
                df.at[row[0], "bonh2"] = row[1]["bonh3"]
                df.at[row[0], "bonh3"] = row[1]["bonh2"]
                df.at[row[0], "bonh4"] = row[1]["bonh1"]
                df.at[row[0], "bonh5"] = row[1]["bonh8"]
                df.at[row[0], "bonh6"] = row[1]["bonh7"]
                df.at[row[0], "bonh7"] = row[1]["bonh6"]
                df.at[row[0], "bonh8"] = row[1]["bonh5"]

                df.at[row[0], "onh1"] = row[1]["onh4"]
                df.at[row[0], "onh2"] = row[1]["onh3"]
                df.at[row[0], "onh3"] = row[1]["onh2"]
                df.at[row[0], "onh4"] = row[1]["onh1"]
                df.at[row[0], "onh5"] = row[1]["onh8"]
                df.at[row[0], "onh6"] = row[1]["onh7"]
                df.at[row[0], "onh7"] = row[1]["onh6"]
                df.at[row[0], "onh8"] = row[1]["onh5"]

                df.at[row[0], "o1"] = row[1]["o4"]
                df.at[row[0], "o2"] = row[1]["o3"]
                df.at[row[0], "o3"] = row[1]["o2"]
                df.at[row[0], "o4"] = row[1]["o1"]
                df.at[row[0], "o5"] = row[1]["o8"]
                df.at[row[0], "o6"] = row[1]["o7"]
                df.at[row[0], "o7"] = row[1]["o6"]
                df.at[row[0], "o8"] = row[1]["o5"]

    df_defrag = gen_feats(df)

    feats = [
        x
        for x in df_defrag.columns.to_list()
        if x
        not in [
            "source",
            "smiles",
            "c_smiles",
            "set",
            "name",
            "ee",
            "ddg",
            "temp",
            "class",
        ]
    ]

    if full:
        return df_defrag
    else:
        return df_defrag[feats]


e_feats = [
    "homo",
    "lumo",
    "nbo1",
    "nbo2",
    "cnbo1",
    "cnbo2",
    "onbo1",
    "onbo2",
    "dipole",
    "nbo-",
    "nbo1/2",
    "gap",
    "nbo*",
    "nbo+",
    "nbom",
    "cnbo-",
    "cnbo1/2",
    "cnbo*",
    "cnbo+",
    "cnbom",
    "onbo-",
    "onbo1/2",
    "onbo*",
    "onbo+",
    "onbom",
]

ta_feats = [
    "na",
    "na_an",
    "nb",
    "estrada",
    "wiener",
    "global_eff",
    "balaban",
    "hosoya",
    "zagreb1",
    "zagreb2",
    "global_simple",
    "kier_a",
    "kier_b",
    "kier_al",
    "0k",
    "1k",
    "2k",
    "3k",
    "1ka",
    "2ka",
    "3ka",
    "1kb",
    "2kb",
    "3kb",
    "1kal",
    "2kal",
    "3kal",
    "k_xia",
    "k_xib",
    "k_xial",
    "redu",
    "0chi",
    "1chi",
    "2chi",
    "3chi",
    "4chi",
    "5chi",
    "0chiv",
    "1chiv",
    "2chiv",
    "3chiv",
    "4chiv",
    "5chiv",
    "T",
    "Tnbo1",
    "Tnbo2",
    "s0chi",
    "s1chi",
    "s2chi",
    "s3chi",
    "s4chi",
    "s5chi",
]

tp_feats = [
    "bpa",
    "crest_flex",
    "angle",
    "t_ova",
    "t_ova_nh",
    "festrada",
    "fglobal_eff",
    "fbalaban",
    "fwiener",
    "fhosoya",
    "fzagreb1",
    "fzagreb2",
    "fkier_a",
    "fkier_b",
    "fkier_al",
    "f0k",
    "f1k",
    "f2k",
    "f3k",
    "f1ka",
    "f2ka",
    "f3ka",
    "f1kb",
    "f2kb",
    "f3kb",
    "f1kal",
    "f2kal",
    "f3kal",
    "k_phi",
    "k_phia",
    "k_phib",
    "k_phial",
    "k_ups",
    "k_upsa",
    "k_upsb",
    "k_upsal",
    "f0chi",
    "f1chi",
    "f2chi",
    "f3chi",
    "f4chi",
    "f5chi",
    "f0chiv",
    "f1chiv",
    "f2chiv",
    "f3chiv",
    "f4chiv",
    "f5chiv",
    "fT",
    "d0chi",
    "d1chi",
    "d2chi",
    "d3chi",
    "d4chi",
    "d5chi",
    "Tenbo1",
    "Tenbo2",
]

sa_feats = [
    "t_vol",
    "t_sur",
    "t_vol_nh",
    "t_sur_nh",
    "b_vol",
    "b_vol_nh",
    "bonh1",
    "bonh2",
    "bonh3",
    "bonh4",
    "bonh5",
    "bonh6",
    "bonh7",
    "bonh8",
    "o1",
    "o2",
    "o3",
    "o4",
    "o5",
    "o6",
    "o7",
    "o8",
    "onh1",
    "onh2",
    "onh3",
    "onh4",
    "onh5",
    "onh6",
    "onh7",
    "onh8",
    "bonh_t",
    "bonhq1",
    "bonhq2",
    "bonhq3",
    "bonhq4",
    "onh_t",
    "onhq1",
    "onhq2",
    "onhq3",
    "onhq4",
    "o_t",
    "oq1",
    "oq2",
    "oq3",
    "oq4",
    "bonhh1",
    "bonhh2",
    "bonhh3",
    "bonhh4",
    "bonhh5",
    "bonhh6",
    "onhh1",
    "onhh2",
    "onhh3",
    "onhh4",
    "onhh5",
    "onhh6",
    "oh1",
    "oh2",
    "oh3",
    "oh4",
    "oh5",
    "oh6",
]

sp_feats = [
    "fbonhq1",
    "fbonhq2",
    "fbonhq3",
    "fbonhq4",
    "fonhq1",
    "fonhq2",
    "fonhq3",
    "fonhq4",
    "foq1",
    "foq2",
    "foq3",
    "foq4",
    "fbonhh1",
    "fbonhh2",
    "fbonhh3",
    "fbonhh4",
    "fbonhh5",
    "fbonhh6",
    "fonhh1",
    "fonhh2",
    "fonhh3",
    "fonhh4",
    "fonhh5",
    "fonhh6",
    "foh1",
    "foh2",
    "foh3",
    "foh4",
    "foh5",
    "foh6",
    "bonhr",
    "onhr",
    "or",
    "pb_vol",
    "pb_vol_nh",
]
