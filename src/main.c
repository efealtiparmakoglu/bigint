/* bigint CLI: fibonacci, factorial and modular exponentiation demos with
 * full decimal output. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bigint.h"

static void usage(void) {
    fprintf(stderr,
            "bigint v1.0.0 — arbitrary-precision arithmetic on base 2^32 limbs\n\n"
            "USAGE\n"
            "  bigint fib <n>               n-th Fibonacci number\n"
            "  bigint fact <n>              n!\n"
            "  bigint powmod <b> <e> <m>    (b ^ e) mod m, decimal arguments\n\n"
            "All values are printed in decimal; arguments may be arbitrarily\n"
            "large (powmod parses them with the same parser).\n\n"
            "Exit codes: 0 ok, 1 failure, 2 usage error\n");
}

int main(int argc, char **argv) {
    if (argc < 3) {
        usage();
        return 2;
    }

    bigint_t x;
    bi_init(&x);

    if (strcmp(argv[1], "fib") == 0) {
        char *end = NULL;
        unsigned long n = strtoul(argv[2], &end, 10);
        if (end == argv[2] || *end != 0) {
            usage();
            return 2;
        }
        bi_fib(&x, (unsigned)n);
    } else if (strcmp(argv[1], "fact") == 0) {
        char *end = NULL;
        unsigned long n = strtoul(argv[2], &end, 10);
        if (end == argv[2] || *end != 0) {
            usage();
            return 2;
        }
        bi_fact(&x, (unsigned)n);
    } else if (strcmp(argv[1], "powmod") == 0 && argc == 5) {
        bigint_t base, exp, mod;
        bi_init(&base);
        bi_init(&exp);
        bi_init(&mod);
        if (bi_from_decimal(&base, argv[2]) != 0 || bi_from_decimal(&exp, argv[3]) != 0 ||
            bi_from_decimal(&mod, argv[4]) != 0) {
            usage();
            return 2;
        }
        bi_mod_pow(&x, &base, &exp, &mod);
        bi_free(&base);
        bi_free(&exp);
        bi_free(&mod);
    } else {
        usage();
        return 2;
    }

    char *text = bi_to_decimal(&x);
    printf("%s\n", text);
    free(text);
    bi_free(&x);
    return 0;
}
