base_folder=/home/ec2-user/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output
dataset_folder=/home/ec2-user/bnet/datasets/user1541667719.91
export NETBOX_HOME=/home/ec2-user/bnetworks_alg/netbox # ../repos/netbox
export PATH=$PATH:$NETBOX_HOME/bin


python $NETBOX_HOME/bin/netAnalyze.py $NETBOX_HOME/conf.props

