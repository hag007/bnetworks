se_folder=/media/hag007/Data/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
dataset_folder=/media/hag007/Data/bnet/datasets/GE_NADAV
export NETBOX_HOME=/home/hag007/repos/bnetworks_alg/netbox # ../repos/netbox
export PATH=$PATH:$NETBOX_HOME/bin

echo $NETBOX_HOME/bin/netAnalyze.py
echo /home/hag007/repos/bnetworks_alg/netbox/conf_0.233509224441.props

python $NETBOX_HOME/bin/netAnalyze.py /home/hag007/repos/bnetworks_alg/netbox/conf_0.233509224441.props


