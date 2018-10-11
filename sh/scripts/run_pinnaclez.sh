base_folder=/home/hag007/bnet
network_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
dataset_folder=$datasets_folder/TNFa
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output


radius=10
max_degree=4000
score=MI
trials=100
improvement=0.001
ST1=0.05
ST2=0.11
ST3=0.00005

java -jar ../repos/pinnaclez/pinnaclez-ORIGINAL.jar $data_folder/TNFa_class \
                                 $data_folder/ge.tsv \
                                 ~/bnet/networks/dip_out.sif \
                                 -r $radius \
                                 -d $max_degree \
                                 -s $score \
                                 -t $trials \
                                 -i $improvement \
                                 -1 $ST1 \
                                 -2 $ST2 \
                                 -3 $ST3 \
                                 -o $output_folder/pinnaclez_results.txt
