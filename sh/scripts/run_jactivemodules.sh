base_folder=/home/hag007/bnet
networks_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
dataset_folder=/home/hag007/bnet/datasets/user1541628615.2
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
is_greedy=True
algo_dir=/home/hag007/repos/bnetworks_alg/jactivemodules

echo /home/hag007/bnet/datasets/user1541628615.2/cache/deg_edger.tsv


java -jar $algo_dir/jactivemodules.jar \
                                 $networks_folder/dip.sif \
                                 /home/hag007/bnet/datasets/user1541628615.2/cache/deg_edger.tsv \
                                 $is_greedy \
                                 /home/hag007/bnet/datasets/user1541628615.2/output/jactivemodules_greedy_results.txt

