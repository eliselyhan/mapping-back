python3 read_star.py -f "../star/run_data_class_1_distinct.star" -c 300.84 310.66 303.87 -e 241.40 349.64 329.70 -x 362.82 285.53 352.82

python3 get_polysomes.py -c run_data_class_1_distinct_tomogram_entry_exit/ -t 150 -o polysome_data/polysomes_run_data_class_1_distinct.csv

python3 write_polysomes_to_star.py -r ../star/run_data_class_1_distinct.star -p polysome_data/polysomes_run_data_class_1_distinct.csv -o class_1_distinct_polysomes_transpose.star
