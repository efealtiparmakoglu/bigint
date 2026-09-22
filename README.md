# bigint

Arbitrary-precision signed integers in **C11** — sign/magnitude over
base-2^32 limbs, with Karatsuba multiplication, Knuth algorithm D
division, modular exponentiation and full decimal conversion. Zero
dependencies, works comfortably at RSA sizes (2048-bit values are part
of the test suite).

```
$ bigint fib 300
222232244629420445529739893461909967206666939096499764990979600...

$ bigint powmod 2 521 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001
4
```

## The pieces

- **Representation** — little-endian `uint32_t` limbs, sign-magnitude
  (zero is always non-negative), capacity-doubling buffers.
- **Multiplication** — schoolbook with 64-bit accumulation, switching to
  **Karatsuba** at 32 limbs (1024 bits): the three recursive products are
  assembled limb-wise so cross-product carries add naturally.
- **Division** — Knuth algorithm D: normalize so the divisor's top limb
  has its high bit set, estimate each quotient digit from the top two
  dividend limbs, correct against the second divisor limb, subtract with
  a documented add-back path. A single-limb fast path handles decimal
  conversion cheaply.
- **modular exponentiation** — square-and-multiply with a reduction after
  every step; scratch buffers never alias the accumulators.
- **Decimal** — parse digit-by-digit (x = x·10 + d); output repeatedly
  divides by 10^9 and prints 9-digit groups. Both directions round-trip
  on 2048-bit values.

## CLI

```
bigint fib <n>               n-th Fibonacci number
bigint fact <n>              n!
bigint powmod <b> <e> <m>    (b ^ e) mod m — decimal args, any size
```

## Tests

```console
$ make test    # 143 checks
```

Highlights: parser round trips up to 2^128-scale strings, sign tables
for add/sub, known 20-digit products, 40 random division identities
(`q·v + r == u`, `|r| < |v|`), a Karatsuba cross-check against an
independent 64-bit modular computation at 2048 bits, Fermat's little
theorem on the Mersenne primes 2^127−1 and 2^521−1, and exact
power-of-two division identities.

## License

MIT — see [LICENSE](LICENSE).
