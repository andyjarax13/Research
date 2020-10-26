"""Empty."""

import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

C = 299792.458  # km/s
DISTANCE = 102.975  # Mpc
H0 = 67.8  # km/s/Mpc
if input("Enter 'c' if conservative. Otherwise realistic.\n") == "c":
    OBJ_SIZE = (0.5, 1, "Conservative")
else:
    OBJ_SIZE = (1, 2, "Realistic")


def create_cluster_data():
    """Empty."""
    cluster_data = []
    cluster_file = open("clusters.csv")
    for line in cluster_file:
        line = line.strip().split(",")
        if not line:
            continue
        if line[0][0] == "#":
            continue
        cluster_data.append(line)
    cluster_file.close()
    cluster_data = pd.DataFrame(
        cluster_data[1:],
        columns=cluster_data[0]).apply(pd.to_numeric, errors="ignore")
    return cluster_data


def physical_separation(ra1, dec1, d1, ra2, dec2, d2):
    """Empty."""
    angle = np.arccos(
        np.sin(dec1 * np.pi / 180) * np.sin(dec2 * np.pi / 180)
        + np.cos(dec1 * np.pi / 180) * np.cos(dec2 * np.pi / 180)
        * np.cos(ra1 * np.pi / 180 - ra2 * np.pi / 180))
    d3 = (d1**2 + d2**2 - 2 * d1 * d2 * np.cos(angle)) ** (1 / 2)
    return d3


def main():
    """Run module."""
    udg_hdul = fits.open("SMUDGes.fit")
    udg_data = pd.DataFrame(udg_hdul[1].data)
    cluster_data = create_cluster_data()

    repeats = {}
    repeats2 = {}
    viable = [",".join(list(cluster_data.columns)) + ",UDG\n"]
    viable2 = [",".join(list(udg_data.columns))
               + ",Distance (mean),Distance (median),Distance (stddev)"
               + ",Size (mean),Size (median),Size (stddev),Cluster\n"]

    for _, udg in udg_data.iterrows():
        d2 = cluster_data["Velocity"] / H0
        d3 = physical_separation(udg["RAJ2000"], udg["DEJ2000"], d2,
                                 cluster_data["RA"], cluster_data["DEC"], d2)
        size = d2 * udg["Re"] * 1e3 / 206265
        filt = ((cluster_data["Type"] == "GClstr") & (d3 < OBJ_SIZE[1]))
        # filt = (((cluster_data["Type"] == "GGroup") & (d3 < OBJ_SIZE[0])) |
        #        ((cluster_data["Type"] == "GClstr") & (d3 < OBJ_SIZE[1])))
        d2 = d2[filt]
        d3 = d3[filt]
        size = size[filt]
        cluster_data_filt = cluster_data[filt]
        for _, cluster in cluster_data_filt.iterrows():
            if str(cluster["Object Name"]) not in repeats:
                viable.append(",".join(map(str, list(cluster))).strip())
                repeats[str(cluster["Object Name"])] = str(udg["SMDG"])
            else:
                repeats[str(cluster["Object Name"])] += "|" + str(udg["SMDG"])
            if str(udg["SMDG"]) not in repeats2:
                viable2.append(
                    ",".join(map(str, list(udg))).strip()
                    + f",{d2.mean()},{d2.median()},{d2.std()}"
                    + f",{size.mean()},{size.median()},{size.std()}")
                repeats2[str(udg["SMDG"])] = str(cluster["Object Name"])
            else:
                repeats2[str(udg["SMDG"])] += "|" + str(cluster["Object Name"])

    for i in range(1, len(viable)):
        key = viable[i].split(",")[0]
        viable[i] += "," + repeats[key] + "\n"
    viable_file = open(f"ViableClusters{OBJ_SIZE[2]}.csv", "w")
    viable_file.writelines(viable)
    viable_file.close()
    for i in range(1, len(viable2)):
        key = viable2[i].split(",")[0]
        viable2[i] += "," + repeats2[key] + "\n"
    viable2_file = open(f"ViableUDGs{OBJ_SIZE[2]}.csv", "w")
    viable2_file.writelines(viable2)
    viable2_file.close()

    fig, ax = plt.subplots(figsize=(10, 8))

    for i in range(1, len(viable)):
        line = viable[i].split(",")
        ra, dec = float(line[1]), float(line[2])
        if line[3] == "GClstr":
            radius = OBJ_SIZE[1] * 2.06265e11 / (
                (float(line[4]) / H0) * 60 * 60 * 1e6)
            col = "g"
        else:
            radius = OBJ_SIZE[0] * 2.06265e11 / (
                (float(line[4]) / H0) * 60 * 60 * 1e6)
            col = "b"
        circle = plt.Circle((ra, dec), radius, color=col, alpha=0.1)
        ax.add_artist(circle)
    ax.plot(udg_data["RAJ2000"], udg_data["DEJ2000"], "r.")
    ax.set_title(f"Sky Around Coma Cluster ({OBJ_SIZE[2]}) "
                 "(Redshift: 0.003 < z < 0.04)")
    ax.set_xlabel("Right Ascension (decimal degrees)")
    ax.set_ylabel("Declination (decimal degrees)")
    ax.invert_xaxis()
    ax.axis("equal")
    ax.legend(loc="lower left", handles=[
        mpatches.Patch(color="r", label="UDG"),
        mpatches.Patch(color="g",
                       label=f"Galaxy Cluster (Radius: {OBJ_SIZE[1]} Mpc)"),
        mpatches.Patch(color="b",
                       label=f"Galaxy Group (Radius: {OBJ_SIZE[0]} Mpc)")
    ])
    fig.tight_layout()
    plt.savefig(f"ViableClusters{OBJ_SIZE[2]}")
    plt.close()


if __name__ == "__main__":
    main()
