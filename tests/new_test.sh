if [[ $# -eq 0 ]] ; then
    echo 'What is the name of the test though?'
    exit 0
fi

mkdir $1
cp trivial/trivial.cpp $1/$1.cpp
sed "s/trivial/$1/g" < trivial/Makemodule.am > $1/Makemodule.am

echo "include tests/$1/Makemodule.am" >> Makemodule.am
