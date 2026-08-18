# Third-party notices

## MaxEntScan model tables

Kea contains a compressed, lossless representation of the MaxEntScan donor and
acceptor lookup tables distributed by the MIT-licensed `maxentpy` project:

- Source: <https://github.com/kepbod/maxentpy>
- Source commit: `b5e69aadd745fcd9ee65c9dd0ff21b902a63a825`
- `score5_matrix.txt` SHA-256:
  `c64fa6aea3d8b71d6af69f5e3ece8f84de06f5450f6d4ba482fbbff6afcf8b5b`
- `score3_matrix.txt` SHA-256:
  `9e8a74dc795ae5c9c5b611d4343a1a5a70fff863956f4e00ae23a85f21beb3ae`
- Bundled `maxent_scan.npz` SHA-256:
  `9f47447284adf66079a51ecb3cadcb2bf7265148e97ed5eba9e0de3f6f9d358f`

The scoring algorithm and tables derive from:

> Yeo G, Burge CB. Maximum entropy modeling of short sequence motifs with
> applications to RNA splicing signals. Journal of Computational Biology.
> 2004;11(2-3):377-394. doi:10.1089/1066527041410418.

The `maxentpy` software and accompanying data are provided under the following
MIT license:

> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.

