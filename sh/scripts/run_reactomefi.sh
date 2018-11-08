base_folder=/home/hag007/bnet
networks_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
dataset_folder=/home/hag007/bnet/datasets/user1541610069.63
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output




java -jar /home/hag007/repos/bnetworks_alg/reactomefi/reactomefi.jar $data_folder/ge.tsv \
                                    $networks_folder/dip.sif \
                                    1.2 \
                                    $output_folder/reactomefi_modules.txt