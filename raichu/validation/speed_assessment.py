from typing import Optional
from scipy.interpolate import make_interp_spline
import numpy as np
from sys import argv
from create_random_clusters import generate_modular_cluster
import time
from statistics import mean, stdev
from dataclasses import dataclass
import matplotlib.pyplot as plt
import os


@dataclass
class SpeedData:
    nr_modules: int
    cluster_type: str
    times: list


def assess_speed(out_folder, max_nr_modules=20):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    module_nrs = []
    nrps_times = []
    nrps_errors = []
    cis_pks_times = []
    cis_pks_errors = []
    trans_pks_times = []
    trans_pks_errors = []
    pks_times = []
    pks_errors = []
    hybrid_times = []
    hybrid_errors = []

    for nr_modules in range(1, max_nr_modules + 1):
        module_nrs.append(nr_modules)
        nrps_speeds = SpeedData(nr_modules, "NRPS", [])
        cis_pks_speeds = SpeedData(nr_modules, "cis-AT PKS", [])
        trans_pks_speeds = SpeedData(nr_modules, "trans-AT PKS", [])
        pks_speeds = SpeedData(nr_modules, "PKS", [])
        hybrid_speeds = SpeedData(nr_modules, "Hybrid", [])

        for i in range(10):
            start_time = time.time()
            generate_modular_cluster(nr_modules, out_folder, nrps=True, cis_pks=False, trans_pks=False,
                                     cluster_nr=i + 100 * nr_modules)
            end_time = time.time()
            time_elapsed = end_time - start_time
            nrps_speeds.times.append(time_elapsed)

            start_time = time.time()
            generate_modular_cluster(nr_modules, out_folder, nrps=False, cis_pks=True, trans_pks=False,
                                     cluster_nr=i + 100 * nr_modules + 10)
            end_time = time.time()
            time_elapsed = end_time - start_time
            cis_pks_speeds.times.append(time_elapsed)

            start_time = time.time()
            generate_modular_cluster(nr_modules, out_folder, nrps=False, cis_pks=False, trans_pks=True,
                                     cluster_nr=i + 100 * nr_modules + 20)
            end_time = time.time()
            time_elapsed = end_time - start_time
            trans_pks_speeds.times.append(time_elapsed)

            start_time = time.time()
            generate_modular_cluster(nr_modules, out_folder, nrps=False, cis_pks=True, trans_pks=True,
                                     cluster_nr=i + 100 * nr_modules + 30)
            end_time = time.time()
            time_elapsed = end_time - start_time
            pks_speeds.times.append(time_elapsed)

            start_time = time.time()
            generate_modular_cluster(nr_modules, out_folder, nrps=True, cis_pks=True, trans_pks=True,
                                     cluster_nr=i + 100 * nr_modules + 40)
            end_time = time.time()
            time_elapsed = end_time - start_time
            hybrid_speeds.times.append(time_elapsed)

        nrps_times.append(mean(nrps_speeds.times))
        nrps_errors.append(stdev(nrps_speeds.times))

        pks_times.append(mean(pks_speeds.times))
        pks_errors.append(stdev(pks_speeds.times))

        cis_pks_times.append(mean(cis_pks_speeds.times))
        cis_pks_errors.append(stdev(cis_pks_speeds.times))

        trans_pks_times.append(mean(trans_pks_speeds.times))
        trans_pks_errors.append(stdev(trans_pks_speeds.times))

        hybrid_times.append(mean(hybrid_speeds.times))
        hybrid_errors.append(stdev(hybrid_speeds.times))

    width = 0.15
    plt.figure(figsize=(10, 6))

    plt.bar([x - 2 * width for x in module_nrs], pks_times, width, yerr=pks_errors, error_kw={"capsize": 2.0, "elinewidth": 1.0}, color='navajowhite', edgecolor="black", label="PKS")
    plt.bar([x - 1 * width for x in module_nrs], cis_pks_times, width, yerr=cis_pks_errors, error_kw={"capsize": 2.0, "elinewidth": 1.0}, color='mistyrose', edgecolor="black", label="cis-AT PKS")
    plt.bar([x + 0 * width for x in module_nrs], trans_pks_times, width, yerr=trans_pks_errors, error_kw={"capsize": 2.0, "elinewidth": 1.0}, color='lemonchiffon', edgecolor="black", label="trans-AT PKS")
    plt.bar([x + 1 * width for x in module_nrs], hybrid_times, width, yerr=hybrid_errors, error_kw={"capsize": 2.0, "elinewidth": 1.0}, color='lightblue', edgecolor="black", label="hybrid")
    plt.bar([x + 2 * width for x in module_nrs], nrps_times, width, yerr=nrps_errors, error_kw={"capsize": 2.0, "elinewidth": 1.0},
            color='mediumseagreen', edgecolor="black", label="NRPS")

    plt.xlabel("# Modules")
    plt.ylabel("Time (s)")
    plt.xticks(module_nrs)

    plt.legend()
    plt.ylim(bottom=0.0)

    # plt.show()
    plt.savefig(os.path.join(out_folder, "speed_assessment.svg"))
    plt.clf()
    plt.figure(figsize=(10, 6))

    smooth_module_nrs = np.linspace(module_nrs[0], module_nrs[-1])

    bspline_pks = make_interp_spline(module_nrs, pks_times)
    smooth_pks = bspline_pks(smooth_module_nrs)

    bspline_cis_pks = make_interp_spline(module_nrs, cis_pks_times)
    smooth_cis_pks = bspline_cis_pks(smooth_module_nrs)

    bspline_trans_pks = make_interp_spline(module_nrs, trans_pks_times)
    smooth_trans_pks = bspline_trans_pks(smooth_module_nrs)

    bspline_nrps = make_interp_spline(module_nrs, nrps_times)
    smooth_nrps = bspline_nrps(smooth_module_nrs)

    bspline_hybrid = make_interp_spline(module_nrs, hybrid_times)
    smooth_hybrid = bspline_hybrid(smooth_module_nrs)

    plt.plot(smooth_module_nrs, smooth_pks, color='navajowhite', label="PKS")
    plt.plot(smooth_module_nrs, smooth_cis_pks, color='mistyrose', label="cis-AT PKS")
    plt.plot(smooth_module_nrs, smooth_trans_pks, color='lemonchiffon', label="trans-AT PKS")
    plt.plot(smooth_module_nrs, smooth_hybrid, color='lightblue', label="hybrid")
    plt.plot(smooth_module_nrs, smooth_nrps, color='mediumseagreen', label="NRPS")

    plt.errorbar(module_nrs, pks_times, yerr=pks_errors,
                 capsize=2.0, elinewidth=1.0, fmt='o', color='navajowhite', ecolor="black")
    plt.errorbar(module_nrs, cis_pks_times, yerr=cis_pks_errors,
                 capsize=2.0, elinewidth=1.0, fmt='o', color='mistyrose', ecolor="black")
    plt.errorbar(module_nrs, trans_pks_times, yerr=trans_pks_errors,
                 capsize=2.0, elinewidth=1.0, fmt='o', color='lemonchiffon', ecolor="black")
    plt.errorbar(module_nrs, hybrid_times, yerr=hybrid_errors,
                 capsize=2.0, elinewidth=1.0, fmt='o', color='lightblue', ecolor="black")
    plt.errorbar(module_nrs, nrps_times, yerr=nrps_errors,
                 capsize=2.0, elinewidth=1.0, fmt='o', color='mediumseagreen', ecolor="black")

    plt.xlabel("# Modules")
    plt.ylabel("Time (s)")
    plt.xticks(module_nrs)

    plt.legend()
    plt.ylim(bottom=0.0)
    plt.savefig(os.path.join(out_folder, "speed_assessment_lineplot.svg"))


if __name__ == "__main__":
    assess_speed(argv[1], int(argv[2]))
