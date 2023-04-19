import re
import warnings


def get_bos(file_name, file_type=None):

    bos = {}

    if file_name[-4:] in [".log", ".out"] or file_type == "gaussian":

        if file_type is None:
            warnings.warn("File assumed to be from Gaussian.")

        tboc_header = r"(?<=\sNatural\sBond\sOrbitals\s\(Summary\):\n\n\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\sPrincipal\sDelocalizations\n\s\s\s\s\s\s\s\s\s\s\sNBO\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\sOccupancy\s\s\s\sEnergy\s\s\s\(geminal,vicinal,remote\)\n\s====================================================================================)"
        tboc_footer = r"(?=\sSorting of NBOs:)"
        nbo_re = r"([0-9 ]{6}\.\s[a-zA-Z* ]{3}.*?)(?=[0-9 ]{6}\.\s[a-zA-Z* ]{3}|\sSorting of NBOs:|-------------------------------)"

        text = open(file_name).read()

        nbocs = re.findall(tboc_header + "(.*?)" + tboc_footer, text, re.DOTALL)
        nbocs_n = [re.findall(nbo_re, x, re.DOTALL) for x in nbocs]

        for nbo in nbocs_n[0]:

            if nbo[8:11] == "BD ":

                if (int(nbo[21:24]) - 1, int(nbo[28:31]) - 1) in bos:
                    bos[int(nbo[21:24]) - 1, int(nbo[28:31]) - 1] += (
                        float(nbo[40:48]) / 2
                    )
                else:
                    bos[int(nbo[21:24]) - 1, int(nbo[28:31]) - 1] = (
                        float(nbo[40:48]) / 2
                    )

            elif nbo[8:11] == "BD*":
                bos[int(nbo[21:24]) - 1, int(nbo[28:31]) - 1] -= float(nbo[40:48]) / 2

                # print(nbo[8:11], int(nbo[12:16]), nbo[18:20], int(nbo[21:24]), nbo[26:28], int(nbo[28:31]), float(nbo[40:48]))

    elif file_name == "wbo" or file_type == "xtb":

        if file_type is None:
            warnings.warn("File assumed to be from XTB.")

        lines = open(file_name).readlines()

        tokens = [line.split() for line in lines]

        bos = {(int(t[0]) - 1, int(t[1]) - 1): float(t[2]) for t in tokens}

    else:
        raise NotImplementedError("File type not recognized or specified.")

    return bos


def get_valences(file_name, file_type=None):

    valences = {}

    if file_name[-4:] in [".log", ".out"] or file_type == "gaussian":

        if file_type is None:
            warnings.warn("File assumed to be from Gaussian.")

        val_header = r"(?<=\sSummary\sof\sNatural\sPopulation\sAnalysis:\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\n\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\n\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\sNatural\sPopulation\s\n\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\sNatural\s\s-----------------------------------------------\n\s\s\s\sAtom\s\sNo\s\s\s\sCharge\s\s\s\s\s\s\s\s\sCore\s\s\s\s\s\sValence\s\s\s\sRydberg\s\s\s\s\s\sTotal\n\s-----------------------------------------------------------------------\n)"
        val_footer = r"(?=\n\s=======================================================================)"

        text = open(file_name).read()

        nvals = re.findall(val_header + "(.*?)" + val_footer, text, re.DOTALL)
        lines = nvals[0].split("\n")

        for line in lines:
            valences[int(line[8:12]) - 1] = float(line[40:48])

    else:
        raise NotImplementedError("File type not recognized or specified.")

    return valences
