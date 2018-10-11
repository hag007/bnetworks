base_folder=/home/hag007/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
dataset_folder=$datasets_folder/TNFa_2
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
is_greedy=true

java -jar ../repos/jactivemodules/jactivemodules.jar \
                                 ~/bnet/networks/dip.sif \
                                 $cache_folder/deg_edger.tsv \
                                 true
