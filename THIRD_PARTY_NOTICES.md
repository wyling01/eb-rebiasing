# Third-Party Notices

## Robust empirical-Bayes confidence intervals (`ebci`)

Parts of the robust EBCI implementation in
`src/pas/estimators/robust_eb.py` and `src/pas/intervals/robust_eb.py`
are Python translations and adaptations of the `R/eb.R` and `R/cv.R`
routines in Michal Kolesár's [`ebci`](https://github.com/kolesarm/ebci)
package.

- Source repository: <https://github.com/kolesarm/ebci>
- Copyright (c) 2021 Michal Kolesár
- License: MIT License

The translated code is used to reproduce the estimator, finite-sample moment
corrections, and robust critical values described by Armstrong, Kolesár, and
Plagborg-Møller (2022), *Robust Empirical Bayes Confidence Intervals*.

## MIT License

Copyright (c) 2021 Michal Kolesár

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
