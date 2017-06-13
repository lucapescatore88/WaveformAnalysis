#! /bin/bash






cd $1

#Remove the header
sed -i '1,24d' C1H*
#Split the data
for f in *.csv; 

do 
	echo "Processing $f file..."; 
	split -d -l 2002 $f $f; 
	
done

rm *.csv
index=0;

for name in *.csv*
do
	cp "${name}" "${index}.csv"
	index=$((index+1))
done


rm -f C1*
