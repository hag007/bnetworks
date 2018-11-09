base_folder=/home/ec2-user/bnet
networks_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
dataset_folder=/home/ec2-user/bnet/datasets/user1541725054.98
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
is_greedy=True
algo_dir=/home/ec2-user/bnetworks_alg/jactivemodules

echo /home/ec2-user/bnet/datasets/user1541725054.98/cache/deg_edger.tsv


java -jar $algo_dir/jactivemodules.jar \
                                 $networks_folder/dip.sif \
                                 /home/ec2-user/bnet/datasets/user1541725054.98/cache/deg_edger.tsv \
                                 $is_greedy \
                                 /home/ec2-user/bnet/datasets/user1541725054.98/output/jactivemodules_greedy_results.txt

