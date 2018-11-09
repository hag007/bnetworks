base_folder=/home/ec2-user/bnet
networks_folder=$base_folder/networks
datasets_folder=$base_folder/datasets
dataset_folder=/home/ec2-user/bnet/datasets/user1541731266.27
data_folder=$dataset_folder/data
cache_folder=$dataset_folder/cache
output_folder=$dataset_folder/output




java -jar /home/ec2-user/bnetworks_alg/reactomefi/reactomefi.jar $data_folder/ge.tsv \
                                    $networks_folder/dip.sif \
                                    1.2 \
                                    $output_folder/reactomefi_modules.txt