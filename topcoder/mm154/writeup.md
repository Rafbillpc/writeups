# Writeup for topcoder MM 154.

[Problem description](https://www.topcoder.com/challenges/7a247ad8-0ae4-47f0-94de-59541be3c469?tab=details)

My solution has two phases. In the first phase (90% of the time), I
try to find a good solution with one large cycle of every color. Due
to the scoring function, most points come from this.

In the second phase (10% of the time), I improve it to an arbitrary
solution.

As a first step, I precompute all intersections between
segments. (This is quite expensive: upto 35% of the computation
time.)

## First phase

The first phase is simulated annealing with the following transitions:
1) Try to insert a new star in one of the cycles. Removing 0/1/2 stars
from a different cycle is allowed.
2) Apply a valid 3-opt move to one of the cycles (not that this does
not change the score, so this is always accepted if valid).

The starting position consists of a random cycle of size 3 for every
color (such that there are no intersections).

Because there is a very large variance / deep local optima, the search
is restarted as many times as possible for 33% of the time. Then the
best 4 states are kept and optimized for 33% of the time. The same
happens for the best 2 states for 33% of the time. (this is similar to
successive halving).

## Second phase

The second phase is hill climbing with a few transitions (the state is
a set of segments):
1) Try to add a segment (removing conflicting segment).
2) Apply a valid 3 opt move to any path/cycle.
3) Try to add all segments in a random order.

## Scores

<details>
<summary>Example scores</summary>
```
01
Score = 405.7277412255662
02
Score = 8998.467120826903
03
Score = 19141.70779646047
04
Score = 9437.63776613454
05
Score = 5874.947430481751
06
Score = 2171.170584987534
07
Score = 8223.765763661806
08
Score = 4808.016459492484
09
Score = 6980.255367821438
10
Score = 12652.08698199629
```
</details>

<details>
<summary>Example scores (time limit x 20)</summary>
```
01
Score = 405.7277412255662
02
Score = 10844.462363734774
03
Score = 20117.347810374973
04
Score = 10249.459117599066
05
Score = 6530.403598620486
06
Score = 2190.3341945320667
07
Score = 9248.77842504568
08
Score = 5187.185668978776
09
Score = 7250.074137551974
10
Score = 12806.464109971963
```
</details>
