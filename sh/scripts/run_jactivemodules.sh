base_folder=/media/hag007/Data/bnet
networks_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
dataset_folder=/media/hag007/Data/bnet/datasets/GE_TNFa_2_recovery_0
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
is_greedy=True
algo_dir=/home/hag007/repos/bnetworks_alg/jactivemodules
num_of_modules=10
overlap_threshold=0

echo /media/hag007/Data/bnet/datasets/GE_TNFa_2_recovery_0/cache/deg_edger.tsv


java -jar $algo_dir/jactivemodules.jar \
                                 /media/hag007/Data/bnet/networks/dip.sif \
                                 /media/hag007/Data/bnet/datasets/GE_TNFa_2_recovery_0/cache/deg_edger.tsv \
                                 $is_greedy \
                                 $num_of_modules \
                                 $overlap_threshold \
                                 /media/hag007/Data/bnet/datasets/GE_TNFa_2_recovery_0/output/jactivemodules_greedy_results.txt

