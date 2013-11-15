regions:
	cc -I/auto/share/include/hdf lcr_regions.c -L/auto/share/lib64 -lmfhdf -ldf -ljpeg -lz -lm -o lcr_regions 

reader:
	gcc -g readDat.c -o readDat

cleanr:
	rm readDat

clean:
	rm lcr_regions test.txt
