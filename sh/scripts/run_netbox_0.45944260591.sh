se_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
dataset_folder=/specific/netapp5/gaga/hagailevi/evaluation/bnet/datasets/GE_random_ERS_1_netbox_256
export NETBOX_HOME=/specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/netbox_0.151908305978 # ../repos/netbox
export PATH=$PATH:$NETBOX_HOME/bin

echo $NETBOX_HOME/bin/netAnalyze.py
echo /specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/netbox_0.151908305978/conf_0.197678989516.props

python $NETBOX_HOME/bin/netAnalyze.py /specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/netbox_0.151908305978/conf_0.197678989516.props


