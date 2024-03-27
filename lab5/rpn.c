#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#define N 10

static void error(unsigned int line, int c, bool *err) {
    if (c >= 32 && c <= 126)
        printf("line %d: error at %c\n", line, (char)c);
    else if (c == '\n')
       printf("line %d: error at \\n\n", line);
    else
        printf("line %d: error at %d\n", line, c);

    *err = true;
}

int main() {
    int stack[N];
    int i = 0, c, x = 0;
    bool num = false, err = false;
    unsigned line = 1;
    int s2;

    while((c = getchar()) != EOF) {

        if (err) {
            if (c == '\n') {
                line++;
                err = false;
                i = x = 0;
            }
            continue;
        }

        if (isdigit(c)) {
            x = x * 10 + c - '0';
            num = true;
            continue;
        }

        if (num) {
            if (i == N) {
                error(line, x, &err);
                continue;
            }
            stack[i++] = x;
            num = false;
            x = 0;
        }

        if (c == '\n') {
            if (i != 1) error(line, c, &err);
            else printf("line %d: %d\n", line, stack[0]);
            line++;
            i = x = 0;
            err = false;
            continue;
        }

        if (strchr("+-*/", c)) {
            if (i < 2) error(line, c, &err);

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

        } else if (!isspace(c) && c != EOF) {
            error(line, c, &err);
        }
    }

    return 0;
}

