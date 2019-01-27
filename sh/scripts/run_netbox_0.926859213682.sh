se_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
dataset_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_TNFa_2_netbox_2418
export NETBOX_HOME=/specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/netbox_0.367283546989 # ../repos/netbox
export PATH=$PATH:$NETBOX_HOME/bin

echo $NETBOX_HOME/bin/netAnalyze.py
echo /specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/netbox_0.367283546989/conf_0.0473448902745.props

python $NETBOX_HOME/bin/netAnalyze.py /specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/netbox_0.367283546989/conf_0.0473448902745.props


