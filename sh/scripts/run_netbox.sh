base_folder=/home/hag007/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
dataset_folder=/home/hag007/bnet/datasets/user1539808512.0
export NETBOX_HOME=/home/hag007/repos/bnetworks_alg/netbox # ../repos/netbox
export PATH=$PATH:$NETBOX_HOME/bin


python $NETBOX_HOME/bin/netAnalyze.py $NETBOX_HOME/conf.props

