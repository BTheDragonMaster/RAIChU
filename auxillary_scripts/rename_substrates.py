import os
from sys import argv


def rename_substrates(cluster_in_file, cluster_out_file,
                      old_substrate_name,
                      new_substrate_name):
    with open(cluster_out_file, 'w') as cluster_out:
        with open(cluster_in_file, "r") as cluster_in:
            header = cluster_in.readline()
            cluster_out.write(header)
            for line in cluster_in:
                line = line.strip()
                if line:
                    line_info = line.split("\t")
                    for i, entry in enumerate(line_info)[:]:
                        if entry == old_substrate_name:
                            line_info[i] = new_substrate_name
                    new_line = '\t'.join(line_info)
                    cluster_out.write(f"{new_line}\n")


def rename_substrates_from_folder(folder, out_folder,
                                  old_substrate_name="METHOXYMALONYL_ACP",
                                  new_substrate_name="METHOXYMALONYL_COA"):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    for cluster_name in os.listdir(folder):
        cluster_folder = os.path.join(folder, cluster_name)
        new_cluster_folder = os.path.join(out_folder, cluster_name)
        if not os.path.exists(new_cluster_folder):
            os.mkdir(new_cluster_folder)
        if os.path.isdir(cluster_folder) and 'cluster' in cluster_name:
            cluster_file_in = os.path.join(cluster_folder, 'cluster.txt')
            cluster_file_out = os.path.join(new_cluster_folder, 'cluster.txt')
            rename_substrates(cluster_file_in, cluster_file_out, old_substrate_name, new_substrate_name)


if __name__ == "__main__":
    rename_substrates_from_folder(argv[1], argv[2])
    