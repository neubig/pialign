pialign
=======

pialign - Phrasal ITG Aligner

`pialign` is a package that allows you to create a phrase table and word alignments from an unaligned parallel corpus.

Author
------

Graham Neubig (neubig at gmail dot com)

Mailing List
------------

If you have a question about pialign, please submit it to the [pialign-users mailing list](http://groups.google.com/group/pialign-users). (If you don't get a response, you can also try contacting Graham at neubig at gmail dot com).

Quick Start
-----------

### Compile

First, `configure` and `make` the program:

    $ autoreconf -i
    $ ./configure
    $ make

#### C++11 Support

This project includes C++11 code and you will need a compliant compiler. You may need to add a specific flag to your compiler, something like -std=c++11, check with your compiler's reference documentation.

### Run

Then, to align a file, say we have a source language file `data/f.txt` and a target language file `data/e.txt`, we can run the program as follows:

    $ mkdir out
    $ pialign data/f.txt data/e.txt out/align. &> out/log.txt &

The program will run for a while (about 1-2 hours for 10k sentences, or 10-20 hours for 100k). Note that by default, only sentences of 40 words or less will be aligned. This can be changed by using `-maxsentlen XX` where XX is the maximum sentence length, but if this is set to a large value the program may use a large amount of memory, and the speed will drop to some extent.

When the program is done running, it will output two files `out/align.1.samp` and `out/align.1.pt`. `align.1.samp` contains the alignments (in the form of trees) created by the aligner, and `align.1.pt` contains a phrase table that can be used with Moses.

### Details/Citation

pialign was written by Graham Neubig while he was an intern at the National Institute of Information and Communication Technology (NICT), Japan.
See `http://www.phontron.com/pialign/` for more documentation.

If you would like more details about the method or want to cite pialign in your research, please reference:

An Unsupervised Model for Joint Phrase Alignment and Extraction
Graham Neubig, Taro Watanabe, Eiichiro Sumita, Shinsuke Mori, and Tatsuya Kawahara
Proceedings of the 49th Annual Meeting of the Association for Computational Linguistics (ACL 2011)


### FAQ

Q: I want to create word alignments.

A: This can be done by running the script `script/itgstats.pl` on `align.1.samp.` There are three types of word alignments that you can create, many-to-many (phrase) alignments, one-to-many (block) alignments, and one-to-one (word) alignments.

  many-to-many:

    $ script/itgstats.pl palign < out/align.1.samp > out/align.1.pal

  one-to-many:

    $ script/itgstats.pl balign < out/align.1.samp > out/align.1.bal

  one-to-one:

    $ script/itgstats.pl align < out/align.1.samp > out/align.1.wal

If you want to visualize these alignments you can do so as follows:

    $ script/visualize.pl data/e.txt data/f.txt out/align.1.pal > out/align.1.vis


Q: I want lexical reordering probabilities for Moses.

A: Run the `itgstats.pl` script on the derivation file. Here, 7 is the maximum length of a phrase, and 0.5 is the `pseudo-count` to be added for smoothing.

    $ script/itgstats.pl lex 7 0.5 < out/align.1.samp | LC_ALL=C sort > out/align.1.lex

  These are bi-directional mono/inv/other probabilities.



Q: What are the scores in the phrase table?

A: These are explained in detail in the referenced paper, but briefly:

- 1,2: The conditional probabilities of the phrases p_t(e|f), p_t(f|e)
- 3: The joint probability of the phrase pair p_t(e,f)
- 4: The average posterior probability of a span containing e,f
- 5,6: Lexical weighting probabilities using model 1 word probabilities
        (only output if the base measure uses model 1)
- 7: The uniform phrase penalty


Q: What is the format of the derivation tree in the .samp file?

A: The bracketed representation in the derivation file has four different types of brackets:

- `[X Y]`: Indicates that children X and Y were generated by a regular ITG symbol
- `<X Y>`: Indicates that children X and Y were generated by a reverse ITG symbol
- `((( e ||| f )))`: Indicates a phrase pair e, f
- `{ X }`: Indicates that X was generated as a single phrase in the model. Words
          inside these brackets were actually aligned together as a single
          phrase, but by default pialign forces alignments down to words to
          preserve word alignments for possible other uses.


Q: pialign is too slow!

A: The easiest way speed up pialign is to use multi-threading. This can be done by setting the
`-batchlen` and `-threads` parameters. `-threads` should be set to the number of cores you want to
use. `-batchlen` must be as large as `-threads`. The larger the value is, the faster multi-threaded
processing will be, but very large values might cause a small decrease in accuracy. It is likely
that a value 10-40 times the number of threads would be appropriate. For example, if you want to use
4 cores, you can set `-threads 4 -batchlen 40`.

Also, if you cannot use multiple cores and you want a speed up at the possible small drop in alignment accuracy, you can reduce the size of the probability beam using the `-probwidth` option. The default is `-probwidth 1e-10` so if you set this to `-probwidth 1e-08` instead, you will likely see a significant speed up.


Q: There is no `configure` file to build the program.

A: If you checked the source directly out from the repository, you may have to run autotools before building. First, make sure autotools is installed on your computer, then run `autoreconf -i` which will prepare the configure file for you.

Options
-------

### Input/Output

Usage:

    $ pialign [OPTIONS] FFILE EFILE PREFIX

- `FFILE` is the foreign input corpus
- `EFILE` is the english input corpus
- `PREFIX` is the prefix that will be used for the output

Other input:

- `-le2f`         A file containing the lexicon probabilities for e2f
- `-lf2e`         A file containing the lexicon probabilities for f2e
               (These can be used with `-base m1` or `-base m1g` but are not necessary)

#### Model Parameters

- `-model`        Model type (hier/len/flat, default: hier)

- `-avgphraselen` A parameter indicating the expected length of a phrase.
               default is small (0.01) to prevent overly long alignments
- `-base`         The type of base measure to use (`m1g` is generally best).
               `m1g`=geometric mean of model 1, `m1`=arithmetic mean of model 1,
               `uni`=simple unigrams (default 'm1g')
- `-defstren`     Fixed strength of the PY process (default none)
- `-defdisc`      Fixed discount of the PY process (default none)
- `-nullprob`    The probability of a null alignment (default 0.01)
- `-noremnull`    Do not remember nulls in the phrase table
- `-termprior`    The prior probability of generating a terminal (0.33)
- `-termstren`    Strength of the type distribution (default 1)

#### Phrase Table

- `-maxphraselen` The maximum length of a minimal phrase (default 7)
- `-maxsentlen`   The maximum length of sentences to use (default 40)
- `-printmax`     The maximum length of phrases included in the phrase table (default 7)
- `-printmin`     The minimal length of phrases included in the phrase table (default 1)
- `-noword`       Output only phrase alignments (do not force output of word alignments)

#### Inference Parameters

- `-burnin`       The number of burn-in iterations (default 9)
- `-histwidth`    The width of the histogram pruning to use (default none)
- `-probwidth`    The width of the probability beam to use (default 1e-10)
- `-samps`        The number of samples to take (default 1)
- `-samprate`     Take samples every samprate turns (default 1)
- `-worditers`    The number of iterations to perform with a word-based model (default 0)
- `-noshuffle`    Don't shuffle the order of the sentences
- `-batchlen`     The number of sentences to process in a single batch
- `-threads`      The number of threads to use (must be <= -batchlen)

Misc.
-----

If you have boost installed, uncomment the lines at the beginning of the `Makefile` to include boost libraries, and uncomment `#define COMPRESS` at the beginning of `definitions.h`. This will make the program automatically compress the output in .gz format.

TODO
----

* Increase the efficiency of Model 1 calculation.
* Reduce the memory used by the phrase table by representing phrases as references to their position in the corpus (like a suffix array).

License
-------

It is distributed under the common public license, version 1.0 (available in the `COPYING` file.)
