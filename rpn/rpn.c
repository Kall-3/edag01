#include <stdio.h>
#include <string.h>
//#include <ctype.h>
//#include <stdbool.h>

#define N 10

static void error(unsigned int line, int c, unsigned int* err) {
    char buf[3];

    if (c == '\n') {
        strcpy(buf, "\\n");
    } else {
        buf[0] = c;
        buf[1] = 0;
    }
    printf("line %u: error at %s\n", line, buf);

    *err = 1;
}

int main() {
    int stack[N];
    int i = 0, c, x = 0;
    unsigned err = 0;
    unsigned line = 1;
    int s2;

    while((c = getchar()) != EOF) {

        if (err) {
            if (c == '\n') {
                line++;
                err = 0;
                i = x = 0;
            }
            continue;
        }

        // if (isdigit(c)) {
        if (c >= 48 && c <= 57) {
            x = x * 10 + c - '0';

            if (i == N) {
                error(line, '0' + x, &err);
                continue;
            }

            while((c = getchar()) >= 48) {
                x = x * 10 + c - '0';
            }
            stack[i++] = x;
            x = 0;
        }

        if (c == '\n') {
            if (i != 1) error(line, c, &err);
            else printf("line %d: %d\n", line, stack[0]);
            line++;
            i = x = 0;
            err = 0;
            continue;
        }

        if (strchr("+-*/", c)) {
            if (i < 2) {
                error(line, c, &err); 
                continue;
            }

            s2 = stack[i-1];

            switch ((char)c) {
                case '+': stack[i-2] += s2; break;
                case '-': stack[i-2] -= s2; break;
                case '*': stack[i-2] *= s2; break;
                case '/': 
                    if (s2 == 0) error(line, c, &err);
                    else stack[i-2] /= s2;
                    break;
            }
            i--;

        }

        if (c == '!') {
            error(line, c, &err);
        }
    }

    return 0;
}

