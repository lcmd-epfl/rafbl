import pandas as pd
from rdkit import Chem


def get_classes(df=None,file=None):

    if file is not None and df is None:
        df = pd.read_csv(file)

    substructs = [
        "C(C/N=C/C)/N=C/C",
        "C1(=NC2C(O1)Cc1c2cccc1)CC1=NC2C(O1)Cc1c2cccc1",
        # "C1(=NC2C(O1)Cc1c2cccc1)NC1=NC2C(O1)Cc1c2cccc1",
        # "C1(=NCCO1)CC1=NCCO1",
        # "C1(=NCCO1)NC1=NCCO1",
        "c1c(nc(cc1)C1=NCCO1)C1=NCCO1",
        "C1(=NCCO1)c1ncccc1",
        "C1(=NCCO1)c1c(cccc1)P",
        "C1CN=C(O1)C1=NCCO1",
        "C1CN=C(N1)C1=NCCN1",
        "C1=NCCO1"
    ]

    classes = []

    for row in df.iterrows():

        mol = Chem.MolFromSmiles(row[1]["c_smiles"])
        class_found = False

        for i,substruct in enumerate(substructs):
            
            if i == 2:
                j = 7
            else:
                j = i
            
            sub = Chem.MolFromSmiles(substruct)
            
            if mol.HasSubstructMatch(sub):
                classes.append(j)
                class_found = True
                break
        
        if not class_found:
            classes.append(-1)


    df["class"] = classes

    return df.copy()
    # df.to_csv("descs2_2.csv",index=False)
    # print(classes,len(classes),len(df.index))
