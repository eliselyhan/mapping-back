#usage: read_star.py [-h] -f FILENAME -c CENTER_MASS CENTER_MASS CENTER_MASS -e
#                    ENTRY ENTRY ENTRY -x EXIT EXIT EXIT -o OUTPUT

#optional arguments:
#  -h, --help            show this help message and exit
#  -f FILENAME, --filename FILENAME
#                        star file with the ribosome data
#  -c CENTER_MASS CENTER_MASS CENTER_MASS, --center_mass CENTER_MASS CENTER_MASS CENTER_MASS
#                        x,y,z coordinates of the subtomogram-averaged
#                        ribosome's center of mass
#  -e ENTRY ENTRY ENTRY, --entry ENTRY ENTRY ENTRY
#                        x,y,z coordinates of the ribosome's entry point
#  -x EXIT EXIT EXIT, --exit EXIT EXIT EXIT
#                        x,y,z coordinates of the ribosome's exit point
#  -o OUTPUT, --output OUTPUT
#                       name for output csv file contatining the entry and
#                        exit coordinates for each ribosomes

python3 read_star.py -f "../star/run_data_class_2.star" -c 300.83 307.31 304.66 -e 252.49 361.03 309.28 -x 357.56 284.28 350.73
