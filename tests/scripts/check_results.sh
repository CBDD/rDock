#!/bin/bash

function get_number_of_atoms {
    # get the number of atoms for the first record of sd file
    # $1 is the file name
    awk 'NR==4 {print $1}' $1
}

function get_coordinates {
    # get the coordinates of the first atom in the file
    # $1 is the file name
    atom_count=$(get_number_of_atoms $1)
    cat $1 | tail -n +5 | head -n ${atom_count} | awk '{print $1, $2, $3}'
}

function get_scores {
    awk '/>  <SCORE>/{getline;print}' $1
}

function compare_scores {
    # compare the final scores of the two files
    # $1 is the first file name
    # $2 is the second file name
    paste <(get_scores $1) <(get_scores $2) | awk '
    function abs(v) {return v < 0 ? -v : v}
        BEGIN {
            # set the tolerance for the comparison
            tolerance = 0.01
        }
        {
            if (abs($1 - $2) > tolerance) {
                printf "[FAIL] scores differ: %f vs %f\n", $1, $2
                exit 1
            } else {
                printf "[OK] scores match\n"
            }
        }
    '
}

function compare_coords {
    # compare the coordinates of the two files
    # $1 is the first file name
    # $2 is the second file name

    paste <(get_coordinates $1) <(get_coordinates $2) | awk '
        function abs(v) {return v < 0 ? -v : v}
        BEGIN {
            # set the tolerance for the comparison
            tolerance = 0.01
            failed = 0
        }
        {
            if (abs($1-$4) > tolerance || abs($2-$5) > tolerance || abs($3-$6) > tolerance) {
                printf "[WARNING] Coordinates differ for atom %i: %f, %f, %f vs %f, %f, %f\n", NR, $1, $2, $3, $4, $5, $6
                failed = 1
            }
        }
        END {
            if (failed) {
                print "[ERROR] coordinates differ"
                exit 1
            } else {
                print "[OK] coordinates match"
                exit 0
            }
        }
    '
}

function compare_results {
    # compare the results of the two files
    # $1 is the first file name
    # $2 is the second file name
    compare_coords $1 $2 || return 1
    compare_scores $1 $2 || return 1
}

function check_results {
    compare_results $1 $2
    if [ $? -eq 0 ]; then
        echo "The test succeeded! The results agree with the reference ones."
        echo "Have fun using rDock!!"
    else
        echo "The test failed, please check the compilation is OK and no errors were raised."
    fi
}

# check this bash script received two arguments. if not, print usage and exit
if [ $# -ne 2 ]; then
    echo "Usage: $0 <sd_file1> <sd_file2>"
    exit 1
fi

compare_results $1 $2
