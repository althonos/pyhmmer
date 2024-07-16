Copyright Notice
================

PyHMMER
-------

The PyHMMER library is developed under the MIT license::

    Copyright (c) 2020-2024 Martin Larralde <martin.larralde@embl.de>

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

PyHMMER development was supported by the European Molecular Biology Laboratory; 
the SFB 1371 of the German Research Foundation (Deutsche Forschungsgemeinschaft, DFG) 
[grant number 395357507] and the Federal Ministry of Education and Research (BMBF) 
[grant number 031L0181A].

Easel
-----

PyHMMER distributes, builds and links to code from the Easel library, 
redistributed under the terms of the BSD license::

    Copyright (C) 1990-2023 Sean R. Eddy
    Copyright (C) 2015-2023 President and Fellows of Harvard College
    Copyright (C) 2000-2023 Howard Hughes Medical Institute
    Copyright (C) 1995-2006 Washington University School of Medicine
    Copyright (C) 1992-1995 MRC Laboratory of Molecular Biology, UK

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
    OF THE POSSIBILITY OF SUCH DAMAGE.

Easel's code includes contributions from members of the Eddy/Rivas
laboratory at Harvard University, the HMMER and Infernal development
teams, and other colleagues, including:

- Tyler Camp
- Nick Carter
- Michael Farrar
- Graeme Mitchison
- Eric Nawrocki
- Sam Petti
- Elena Rivas
- Travis Wheeler

Easel also includes code we have incorporated from other sources and
authors -- including public domain code, and licensed copyrighted
code. Sources and licenses are noted in the appropriate places in
individual files. Copyright holders and contributors include:

- Barry W. Brown, James Lovato:        ``esl_random:esl_rnd_Gaussian()``
- Bob Jenkins:                         ``esl_random::esl_rnd_mix3()``
- Steven G. Johnson, and others:       autoconf macros in m4/
- Martin Larralde:                     ARM Neon support
- Kevin Lawler:                        ``esl_rand64::esl_rand64_Deal()``
- Stephen Moshier:                     SIMD vectorized ``logf``, ``expf``
- Takuji Nishimura, Makoto Matsumoto:  ``esl_random``, ``esl_rand64``
- Julien Pommier:                      SIMD vectorized ``logf``,``expf``
- David Robert Nadeau:                 ``esl_stopwatch``
- Henry Spencer:                       ``esl_regexp``
- David Wheeler:                       ``easel::esl_tmpfile()``
- Free Software Foundation, Inc.:      ``configure``
- FreeBSD:                             ``easel::esl_strsep()``
- Sun Microsystems, Inc.:              ``esl_stats::esl_erfc()``

Easel development is supported in part by the National Human Genome
Research Institute of the US National Institutes of Health under grant
number R01HG009116. The content is solely the responsibility of the
authors and does not necessarily represent the official views of the
National Institutes of Health.


HMMER
-----

PyHMMER distributes, builds and links to code from the HMMER software,
redistributed under the terms of the BSD-3-clause license::

    Copyright (C) 1992-2023 Sean R. Eddy
    Copyright (C) 2015-2023 President and Fellows of Harvard College
    Copyright (C) 2000-2023 Howard Hughes Medical Institute
    Copyright (C) 1995-2006 Washington University School of Medicine
    Copyright (C) 1992-1995 MRC Laboratory of Molecular Biology

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions
    are met:

    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

    3. Neither the name of any copyright holder nor the names of
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
    OF THE POSSIBILITY OF SUCH DAMAGE.

The code includes contributions and input from current and past
members of the HMMER development team, as well as other colleagues and
sources, including:

- Bill Arndt
- Jeremy Buhler
- Tyler Camp
- Nick Carter
- Sergi Castellano
- Goran Ceric
- Michael Farrar
- Rob Finn
- Ian Holmes
- Bjarne Knudsen
- Diana Kolbe
- Martin Larralde
- Erik Lindahl
- Graeme Mitchison
- Eric Nawrocki
- Lee Newberg
- Sam Petti
- Elena Rivas
- Walt Shands
- Travis Wheeler

HMMER also includes copyrighted and licensed code that has been
incorporated from other sources, including:

- Yuta Mori (libdivsufsort-lite)
- Apple Computer
- Free Software Foundation, Inc.
- IBM TJ Watson Research Center
- X Consortium

HMMER uses the Easel software library, which has its own license and
copyright information. See `Easel`_ section above.

HMMER includes patent-pending SIMD technology under a nonexclusive
license from the estate of Michael Farrar. You are sublicensed to use
this technology specifically for the use, modification, and
redistribution of HMMER.

HMMER development is supported in part by the National Human Genome
Research Institute of the US National Institutes of Health under grant
number R01HG009116. The content is solely the responsibility of the
authors and does not necessarily represent the official views of the
National Institutes of Health.


