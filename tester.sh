# bash file to automatically run all tests in the test suite folder
for in_f in tests/*.in
do
    out_f="$(echo $in_f | cut -d '.' -f 1).out"
    ./pagerank < $in_f | diff - $out_f
    if [ $? -eq "1" ]
    then
        echo "The test $in_f has failed!"
    else
        echo "The test $in_f has passed!"
    fi
done
