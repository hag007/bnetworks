base_folder=/home/ec2-user/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
dataset_folder=/home/ec2-user/bnet/datasets/user1541677281.14

pwd

java -jar -Xmx2G KPM-4.0.jar -strategy=INES -algo=GREEDY -K=1  -L1=420
