se_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
dataset_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_ROR_1_netbox_443
export NETBOX_HOME=/specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/netbox_0.264498371273 # ../repos/netbox
export PATH=$PATH:$NETBOX_HOME/bin

echo $NETBOX_HOME/bin/netAnalyze.py
echo /specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/netbox_0.264498371273/conf_0.650943568799.props

python $NETBOX_HOME/bin/netAnalyze.py /specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/netbox_0.264498371273/conf_0.650943568799.props


