from raichu.run_raichu import read_in_cluster, draw_cluster, draw_product
import os

dirname_MIBiG_BGC0000178 = "MIBiG/Trans_AT_PKS/Cluster_outputs/BGC0000178"

cluster_BGC0000178, cluster_representation_BGC0000178 = read_in_cluster(
    dirname_MIBiG_BGC0000178
)

draw_cluster(
    cluster_representation_BGC0000178,
    out_file=os.path.join(dirname_MIBiG_BGC0000178, "cluster.svg"),
)
draw_product(
    cluster_representation_BGC0000178,
    out_file=os.path.join(dirname_MIBiG_BGC0000178, "final_product.svg"),
)
