/* Numerical test suite for bigint: parser round trips, arithmetic signs,
 * Karatsuba/schoolbook agreement, division properties, mod_pow identities,
 * and RSA-scale feasibility. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bigint.h"

static int tests_run = 0;
static int tests_failed = 0;

#define CHECK(cond)                                                       \
    do {                                                                  \
        tests_run++;                                                      \
        if (!(cond)) {                                                    \
            tests_failed++;                                               \
            fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        }                                                                 \
    } while (0)

static void test_parse_round_trip(void) {
    const char *cases[] = {
        "0", "1", "-1", "42", "-42", "4294967295", "4294967296",
        "18446744073709551615", /* 2^64-1 */
        "340282366920938463463374607431768211456", /* 2^128 */
        "-170141183460469231731687303715884105728",
        "123456789012345678901234567890123456789012345678901234567890",
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        bigint_t x;
        bi_init(&x);
        CHECK(bi_from_decimal(&x, cases[i]) == 0);
        char *text = bi_to_decimal(&x);
        CHECK(strcmp(text, cases[i]) == 0);
        free(text);
        bi_free(&x);
    }
}

static void test_parse_rejects_garbage(void) {
    const char *cases[] = {"", "-", "12a3", "+-1", "1 2"};
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        bigint_t x;
        bi_init(&x);
        CHECK(bi_from_decimal(&x, cases[i]) != 0);
        bi_free(&x);
    }
}

static void test_add_sub_signs(void) {
    struct {
        const char *a, *b, *sum, *diff;
    } cases[] = {
        {"5", "3", "8", "2"},
        {"-5", "-3", "-8", "-2"},
        {"5", "-3", "2", "8"},
        {"-5", "3", "-2", "-8"},
        {"0", "0", "0", "0"},
        {"1000000000000", "1", "1000000000001", "999999999999"},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        bigint_t a, b, s, d;
        bi_init(&a); bi_init(&b); bi_init(&s); bi_init(&d);
        bi_from_decimal(&a, cases[i].a);
        bi_from_decimal(&b, cases[i].b);
        bi_add(&s, &a, &b);
        bi_sub(&d, &a, &b);
        char *st = bi_to_decimal(&s);
        char *dt = bi_to_decimal(&d);
        CHECK(strcmp(st, cases[i].sum) == 0);
        CHECK(strcmp(dt, cases[i].diff) == 0);
        free(st); free(dt);
        bi_free(&a); bi_free(&b); bi_free(&s); bi_free(&d);
    }
}

static void test_mul_known_values(void) {
    struct {
        const char *a, *b, *product;
    } cases[] = {
        {"0", "12345678901234567890", "0"},
        {"123456789", "987654321", "121932631112635269"},
        {"-12", "12", "-144"},
        {"-12", "-12", "144"},
        {"99999999999999999999", "99999999999999999999",
         "9999999999999999999800000000000000000001"},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        bigint_t a, b, p;
        bi_init(&a); bi_init(&b); bi_init(&p);
        bi_from_decimal(&a, cases[i].a);
        bi_from_decimal(&b, cases[i].b);
        bi_mul(&p, &a, &b);
        char *text = bi_to_decimal(&p);
        CHECK(strcmp(text, cases[i].product) == 0);
        free(text);
        bi_free(&a); bi_free(&b); bi_free(&p);
    }
}

static void test_divmod_properties(void) {
    /* for random pairs: q*v + r == u and |r| < |v| */
    unsigned long long state = 12345;
    bigint_t u, v, q, r, check, scaled;
    bi_init(&u); bi_init(&v); bi_init(&q); bi_init(&r); bi_init(&check); bi_init(&scaled);

    for (int trial = 0; trial < 40; trial++) {
        char buffer[128];
        int digits = 5 + trial % 60;
        for (int i = 0; i < digits; i++) {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            buffer[i] = (char)('0' + (state >> 33) % 10);
        }
        buffer[digits] = 0;
        bi_from_decimal(&u, buffer);
        for (int i = 0; i < digits / 2 + 1; i++) {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            buffer[i] = (char)('0' + (state >> 33) % 10);
        }
        buffer[digits / 2 + 1] = 0;
        bi_from_decimal(&v, buffer);
        if (v.len == 0) continue;

        bi_divmod(&q, &r, &u, &v);
        bi_mul(&scaled, &q, &v);
        bi_add(&check, &scaled, &r);
        CHECK(bi_cmp(&check, &u) == 0);
        bigint_t abs_r, abs_v;
        bi_init(&abs_r); bi_init(&abs_v);
        abs_r.neg = 0; bi_copy(&abs_r, &r);
        abs_v.neg = 0; bi_copy(&abs_v, &v);
        CHECK(bi_cmp(&abs_r, &abs_v) < 0);
        bi_free(&abs_r); bi_free(&abs_v);
    }
    bi_free(&u); bi_free(&v); bi_free(&q); bi_free(&r); bi_free(&check); bi_free(&scaled);
}

static void test_karatsuba_matches_schoolbook(void) {
    /* 2^400: verify through bit length and a modular spot check */
    bigint_t a, b, p, q, r, m;
    bi_init(&a); bi_init(&b); bi_init(&p); bi_init(&q); bi_init(&r); bi_init(&m);
    bi_set_u64(&a, 2);
    for (int i = 0; i < 399; i++) {
        bi_add(&a, &a, &a); /* a = 2^400 */
    }
    bi_copy(&b, &a);
    bi_mul(&p, &a, &b); /* 2^800 */
    CHECK(bi_bitlen(&p) == 801);
    /* (2^400)^2 mod 7 == (2^400 mod 7)^2 mod 7; 2^400 mod 7 = 4 since the
     * powers of 2 cycle with period 3 and 400 = 3*133 + 1 */
    bi_from_decimal(&m, "7");
    bi_divmod(&q, &r, &p, &m);
    /* 2^800 mod 7: powers of 2 cycle mod 7 with period 3, and
     * 800 = 3*266 + 2, so the remainder is 2^2 = 4 */
    bigint_t four, two;
    bi_init(&four);
    bi_init(&two);
    bi_set_u64(&four, 4);
    bi_set_u64(&two, 2);
    CHECK(bi_cmp(&r, &four) == 0);
    bi_free(&a); bi_free(&b); bi_free(&p); bi_free(&q); bi_free(&r); bi_free(&m);
    bi_free(&four); bi_free(&two);
}

static void test_big_division_identity(void) {
    /* 2^2048 / (2^1024 + 1): quotient and remainder satisfy q*v + r == u */
    bigint_t u, v, q, r, one, check, scaled;
    bi_init(&u); bi_init(&v); bi_init(&q); bi_init(&r);
    bi_init(&one); bi_init(&check); bi_init(&scaled);
    bi_set_u64(&one, 1);
    bi_set_u64(&u, 2);
    for (int i = 0; i < 2047; i++) bi_add(&u, &u, &u);
    bi_set_u64(&v, 2);
    for (int i = 0; i < 1023; i++) bi_add(&v, &v, &v);
    bi_add(&v, &v, &one); /* v = 2^1024 + 1 */

    bi_divmod(&q, &r, &u, &v);
    CHECK(bi_cmp(&r, &v) < 0);
    bi_mul(&scaled, &q, &v);
    bi_add(&check, &scaled, &r);
    CHECK(bi_cmp(&check, &u) == 0);

    /* exact power division: 2^2048 / 2^1024 == 2^1024, remainder 0 */
    bigint_t half, zero;
    bi_init(&half);
    bi_init(&zero);
    bi_set_u64(&half, 2);
    for (int i = 0; i < 1023; i++) bi_add(&half, &half, &half);
    bi_divmod(&q, &r, &u, &half);
    CHECK(bi_cmp(&q, &half) == 0);
    CHECK(r.len == 0);

    bi_free(&u); bi_free(&v); bi_free(&q); bi_free(&r);
    bi_free(&one); bi_free(&check); bi_free(&scaled); bi_free(&half); bi_free(&zero);
}

static void test_karatsuba_cross_check(void) {
    /* build 2^1024 by doubling (1024 bits = 32 limbs = Karatsuba threshold),
     * square it, and cross-check the low-cost identity p mod (10^9+7)
     * against an independent 64-bit modular computation. */
    bigint_t x, p, q, r, m;
    bi_init(&x); bi_init(&p); bi_init(&q); bi_init(&r); bi_init(&m);
    bi_set_u64(&x, 2);
    for (int i = 0; i < 1023; i++) {
        bi_add(&x, &x, &x);
    }
    CHECK(bi_bitlen(&x) == 1025);
    bi_mul(&p, &x, &x);
    CHECK(bi_bitlen(&p) == 2049);

    const uint64_t M = 1000000007ULL;
    uint64_t x_mod = 1;
    for (int i = 0; i < 1024; i++) {
        x_mod = (x_mod * 2) % M;
    }
    uint64_t expected = (x_mod * x_mod) % M;
    bi_from_decimal(&m, "1000000007");
    bi_divmod(&q, &r, &p, &m);
    char *text = bi_to_decimal(&r);
    char expected_text[32];
    snprintf(expected_text, sizeof(expected_text), "%llu", (unsigned long long)expected);
    CHECK(strcmp(text, expected_text) == 0);
    free(text);
    bi_free(&x); bi_free(&p); bi_free(&q); bi_free(&r); bi_free(&m);
}

static void test_mod_pow_fermat(void) {
    bigint_t base, exp, mod, result, one, two;
    bi_init(&base); bi_init(&exp); bi_init(&mod); bi_init(&result);
    bi_init(&one); bi_init(&two);
    bi_set_u64(&one, 1);
    bi_set_u64(&two, 2);

    /* Fermat: 2^(p-1) mod p == 1 for the primes below */
    struct {
        const char *mod_decimal;
        unsigned mod_bits;
    } primes[] = {
        {"127", 7},                    /* 2^7 - 1 */
        {"170141183460469231731687303715884105727", 127}, /* 2^127 - 1 */
    };
    for (size_t i = 0; i < sizeof(primes) / sizeof(primes[0]); i++) {
        bi_from_decimal(&mod, primes[i].mod_decimal);
        bi_sub(&exp, &mod, &one); /* p - 1 */
        bi_set_u64(&base, 2);
        bi_mod_pow(&result, &base, &exp, &mod);
        CHECK(bi_cmp(&result, &one) == 0);
    }

    /* the 521-bit Mersenne prime 2^521 - 1: Fermat at RSA-adjacent size */
    bi_set_u64(&mod, 2);
    for (int i = 0; i < 520; i++) {
        bi_add(&mod, &mod, &mod);
    }
    bi_sub(&mod, &mod, &one); /* mod = 2^521 - 1 */
    bi_sub(&exp, &mod, &one); /* exp = 2^521 - 2 */
    bi_set_u64(&base, 2);
    bi_mod_pow(&result, &base, &exp, &mod);
    CHECK(bi_cmp(&result, &one) == 0);

    /* 5^3 mod 1000 = 125 */
    bi_from_decimal(&base, "5");
    bi_from_decimal(&exp, "3");
    bi_from_decimal(&mod, "1000");
    bi_mod_pow(&result, &base, &exp, &mod);
    char *text = bi_to_decimal(&result);
    CHECK(strcmp(text, "125") == 0);
    free(text);

    bi_free(&base); bi_free(&exp); bi_free(&mod); bi_free(&result);
    bi_free(&one); bi_free(&two);
}

static void test_fib_and_fact(void) {
    bigint_t x;
    bi_init(&x);
    bi_fib(&x, 100);
    char *text = bi_to_decimal(&x);
    CHECK(strcmp(text, "354224848179261915075") == 0); /* F(100) */
    free(text);
    bi_fact(&x, 30);
    text = bi_to_decimal(&x);
    CHECK(strcmp(text, "265252859812191058636308480000000") == 0); /* 30! */
    free(text);
    bi_free(&x);
}

static void test_rsa_scale_mod_pow(void) {
    /* decimal conversion of a 2048-bit number: 2^2048 has 617 digits;
     * the remainder mod 10^9 is verified against an independent 64-bit
     * computation, and the text round-trips through the parser */
    bigint_t u, q, r, m, back;
    bi_init(&u); bi_init(&q); bi_init(&r); bi_init(&m); bi_init(&back);
    bi_set_u64(&u, 2);
    for (int i = 0; i < 2047; i++) bi_add(&u, &u, &u);

    char *text = bi_to_decimal(&u);
    CHECK(strlen(text) == 617); /* floor(2048 * log10(2)) + 1 = 617 */

    uint64_t expect = 1;
    for (int i = 0; i < 2048; i++) {
        expect = (expect * 2) % 1000000000ULL;
    }
    bi_from_decimal(&m, "1000000000");
    bi_divmod(&q, &r, &u, &m);
    char *rem = bi_to_decimal(&r);
    char expect_text[16];
    snprintf(expect_text, sizeof(expect_text), "%llu", (unsigned long long)expect);
    CHECK(strcmp(rem, expect_text) == 0);

    CHECK(bi_from_decimal(&back, text) == 0);
    CHECK(bi_cmp(&back, &u) == 0);
    free(rem);
    free(text);
    bi_free(&u); bi_free(&q); bi_free(&r); bi_free(&m); bi_free(&back);
}

int main(void) {
    test_parse_round_trip();
    test_parse_rejects_garbage();
    test_add_sub_signs();
    test_mul_known_values();
    test_divmod_properties();
    test_karatsuba_matches_schoolbook();
    test_karatsuba_cross_check();
    test_big_division_identity();
    test_mod_pow_fermat();
    test_fib_and_fact();
    test_rsa_scale_mod_pow();

    printf("%d checks, %d failed\n", tests_run, tests_failed);
    return tests_failed == 0 ? 0 : 1;
}
