base_folder=/home/ec2-user/bnet
networks_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
<<<<<<< HEAD
dataset_folder=/home/ec2-user/bnet/datasets/user1541725054.98
=======
dataset_folder=/home/hag007/bnet/datasets/user1541722326.25
>>>>>>> 153e767b3d5b5fd24e3da9f61a25665639360a27
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
is_greedy=True
algo_dir=/home/ec2-user/bnetworks_alg/jactivemodules

<<<<<<< HEAD
echo /home/ec2-user/bnet/datasets/user1541725054.98/cache/deg_edger.tsv
=======
echo /home/hag007/bnet/datasets/user1541722326.25/cache/deg_edger.tsv
>>>>>>> 153e767b3d5b5fd24e3da9f61a25665639360a27


java -jar $algo_dir/jactivemodules.jar \
                                 $networks_folder/dip.sif \
<<<<<<< HEAD
                                 /home/ec2-user/bnet/datasets/user1541725054.98/cache/deg_edger.tsv \
                                 $is_greedy \
                                 /home/ec2-user/bnet/datasets/user1541725054.98/output/jactivemodules_greedy_results.txt
=======
                                 /home/hag007/bnet/datasets/user1541722326.25/cache/deg_edger.tsv \
                                 $is_greedy \
                                 /home/hag007/bnet/datasets/user1541722326.25/output/jactivemodules_greedy_results.txt
>>>>>>> 153e767b3d5b5fd24e3da9f61a25665639360a27

