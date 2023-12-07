#include <stdio.h>
#include <stdlib.h>

typedef struct poly_t {
    int a;
    int b;
} poly_t;

int main(int argc, char *argv[])
{
    char* arr;
    arr = calloc(100, sizeof(char));
    arr[0] = 7;
    arr[1] = 'a';

    printf("%s\n", arr);
}
