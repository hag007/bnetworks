base_folder=/home/hag007/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
dataset_folder=$datasets_folder/TNFa_2
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
is_greedy=true

pwd

java -jar -Xmx2G KPM-4.0.jar -strategy=INES -algo=GREEDY -K=1  -L1=420
