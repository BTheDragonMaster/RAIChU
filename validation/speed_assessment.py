from create_random_clusters import generate_modular_cluster
import time
from statistics import mean


def assess_speed(nr_modules=10, nr_clusters=1000, nrps=True, cis_pks=True, trans_pks=True):
    all_times = []

    for i in range(1000):
        start_time = time.time()
        generate_modular_cluster(nr_modules, nrps=nrps, cis_pks=cis_pks, cluster_nr=i)
        end_time = time.time()
        time_elapsed = end_time - start_time
        all_times.append(time_elapsed)

    return mean(all_times)

if __name__ == "__main__":
