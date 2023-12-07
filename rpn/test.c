#include <stdio.h>

int main(int argc, char *argv[])
{
    int stack[10];
    int i = 0;

    stack[0] = 1;
    stack[1] = 2;

    int s1 = stack[i];
    int s2 = stack[i+1];

    s1 += s2;

    printf("%d\n", s1);

    return 0;
}
