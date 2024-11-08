# link2Rhtslib
To ease testing of proper linkage to the HTSlib library in [`Rhtslib`](https://github.com/nhayden/Rhtslib).

- `callLocal`: native routine with no dependencies
- `htsVersion`: native routine that requires `hts.h` header from HTSlib and links to `libhts` in [`Rhtslib`](https://github.com/nhayden/Rhtslib)
