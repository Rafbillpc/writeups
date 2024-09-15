# Writeup for AHC037



[Problem description](https://atcoder.jp/contests/ahc037/tasks/ahc037_a)
[Standings](https://atcoder.jp/contests/ahc037/standings)
[Best submission](TODO)

The solution can be encoded as a binary tree, with the input points as leaves.
The inner nodes correspond to intermediate points that are constructed to reach the input points more efficiently.
Given the shape of the binary tree, we can compute the positions of the inner nodes: if its children are located at $$(x_1,y_1)$$ and $$(x_2,y_2)$$,
then the inner node is located at $$(min(x_1,x_2),min(y_1,y_2))$$.

The total cost of the binary tree is then the sum of $$\|x_1-x_2\|+\|y_1-y_2\|$$ over all inner nodes. This value is the cost of the inner node.

This total cost is minimized using simulated annealing (the initialization does not matter much). 
There are two transitions:

1. (~85% of all transitions) Move a node `x₁` as a new neighbor of another node `x₂`.
```
      q₁           q₂                                 q₂
     /            /                     q₁           /
    p₁           p₂          ===>      /            p₂
   /  \         /  \                  y₁           /  \
  x₁   y₁      x₂   y₂                            p₁   y₂
                                                 / \
                                                x₁  x₂
```
This is only possible if `x₁` is not an ancestor of `x₂`, and `x₂` is not an ancestor of `x₁`.

2. (~15% of all transitions) Do a tree rotation.
```
        r₁                    r₁
       /                       \
      q₁                        p₁
     / \         ===>          / \ 
    p₁  p₂                    x₁  q₁
   /  \                          /  \
  x₁   y₂                       x₂   p₂
```

These transitions can be represented as a sequence of edge deletions and edge creations.
After performing these changes, we need to recompute the positions and costs of some of the inner nodes.

We can simply move up the tree and recompute this data, until we encounter an inner node whose position stays the same.

With a careful implementation, up to ~50-60 millions of iterations can be done in 2 seconds.

It was quite important to prune bad moves for the first transition: if the nodes x₁ and x₂ are too far away, we don't evaluate the transition.
