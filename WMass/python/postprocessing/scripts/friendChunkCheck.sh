#!/bin/bash
if [[ "$1" == "" || "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 [ -z ] [ dir ]"
    echo "Will look for chunks and check if there is some one missing. "
    echo "It only looks for holes in the sequence of chunks, " 
    echo "and the presence of zombies if option -z is given. " 
    exit 1;
fi

Z=0
if [[ "$1" == "-z" ]]; then
    echo "# Will also check if rootfiles $F are zombies or not"
    Z=1; shift;
fi;

dir="./"
if [[ "$1" != "" ]]; then dir=$1; fi
for F in $(ls ${dir}/*_Friend_*.chunk*.root | sed 's/\.chunk[0-9]\+//' | sort | uniq); do
    FILES=$(ls ${F/.root/.chunk*.root} | \
            perl -npe 's/\.chunk(\d+)\./sprintf(".%06d.",$1)/e' | \
            sort -n | \
            perl -npe 's/\.(\d+)\.root$/sprintf(".chunk%d.root",$1)/e' );
    # echo -e "\nCheck chunk files for $F"; 
    filesimple=$(ls ${F/.root/.chunk*.root})
    filesarray=($FILES)
    NCHUNKS=$((${#filesarray[@]} - 1 ))
    for c in `seq 0 $NCHUNKS`; do
        ftest=$(echo $F | awk -F "." '{print $1 ".chunk"}');
        ftest2=$ftest$c".root"
        if [ ! -f $ftest2 ]; then
            echo "$ftest2 # not present";
        else
            ftest3=$(stat --printf="%s" $ftest2)
            if [ $ftest3 -lt 1000 ]; then
                echo "$ftest2 # has zero size";
            fi;
        fi;
    done
done

if [[ "$Z" != "0" ]]; then
    echo "# Testing for zombies or not correctly closed files";
    FILES=$(ls ${dir}/*_Friend_*.chunk*.root);
    id=`shuf -i 1-100 -n 1`
    for Z in $FILES; do
        if test -s $Z; then # empty files have already been found
            root -b -l -q $Z >& check_$id.log 
            result=$(grep -E "(nullptr|recover|Zombie)" check_$id.log | wc -l)
            if [ $result -ne 0 ]; then
                echo "$Z     # zombie";
            else
                echo "$Z     # OK";
            fi;
            rm check_$id.log
        fi;
    done
fi;

