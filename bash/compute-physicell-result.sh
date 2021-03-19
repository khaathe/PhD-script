#/usr/bin/bash


result_dir="./result"
log_file="$result_dir/log.txt"

if [ $# -eq 0 ]; then
    last_povwriter_index=`ls output/ | grep -c "output0.*.xml"`
    let "last_povwriter_index=$last_povwriter_index-1"
else
    last_povwriter_index=$1
fi

if [ -d "$result_dir/pov" ]; then 
    rm $result_dir/pov/*.pov
fi

if [ -d "$result_dir/png" ]; then 
    rm $result_dir/png/*.png
fi

if [ -d "$result_dir/gif" ]; then 
    rm $result_dir/gif/*.gif
fi

echo ">>>>> mkdir $result_dir"
echo ">>>>> mkdir $result_dir/pov $result_dir/png $result_dir/gif"

echo ">>>>> povwriter 0:1:$last_povwriter_index"
povwriter 0:1:$last_povwriter_index

echo ">>> mv *.pov $result_dir/pov"
mv *.pov "$result_dir/pov"

for  pov_file in `ls $result_dir/pov`; do
    png_file=$(echo "$pov_file"|sed 's/pov//g')
    echo ">>> creating png : $result_dir/png/$png_file"
    povray +O$result_dir/png/$png_file  $result_dir/pov/$pov_file
done

echo "magick convert $result_dir/png/*.png $result_dir/gif/simu.gif"
magick convert $result_dir/png/*.png $result_dir/gif/simu.gif
