directory=$1

for f in $directory*; 
do
  b=$(basename $f)
  #echo $f $b
  ./FLASH_Analysis.bash $f/ $b & 
done
wait
