for x in s 0 1 2 3 4 5
do
	gcc -c -O$x -o $x.o mipsx.c
done
