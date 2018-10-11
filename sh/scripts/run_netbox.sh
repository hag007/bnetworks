base_folder=/home/hag007/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
dataset_folder=$datasets_folder/TNFa_2
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
export NETBOX_HOME=../repos/netbox
export PATH=$PATH:$NETBOX_HOME/bin


python $NETBOX_HOME/bin/netAnalyze.py $cache_folder/netbox1.props