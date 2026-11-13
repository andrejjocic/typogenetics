# typogenetics

An implementation of *typogenetics*, the toy model of molecular biology from Douglas Hofstadter's *Gödel, Escher, Bach* (GEB), together with a program that searches for **self-replicating strands**: the puzzle Hofstadter left to the reader.

This README summarises my bachelor's thesis, [Searching for self-replicators in typogenetics](https://repozitorij.uni-lj.si/IzpisGradiva.php?id=139590&lang=eng) (written in Slovenian).

## Usage

Requires **Python 3.10+** (structural pattern matching). `search.py` also needs `matplotlib` and `numpy`.

```sh
python interpreter.py                    # execute typo-acids one at a time on a strand of your choice
python search.py -c 20 hrs -p 3          # heuristic search (see Results)
python search.py -c 20 hrs -p 2 -s 15    # heuristic search with the similarity-based fitness
python search.py -c 13 lex               # lexicographic search
```

The search runs until you press Ctrl+C, then prints all self-replicators found. Useful options:

| option | meaning |
|---|---|
| `-c L` | only look for cycles of length ≤ L in the daughter graph (strongly recommended) |
| `-v`, `-vv`, `-vvv` | report each new batch of self-replicators / also print the cycles proving it / also log every expanded strand |
| `-p` (before `lex`/`hrs`) | plot search statistics when done |
| `hrs -p P` | expand P queued strands for every fresh one (default 1) |
| `hrs -s W` | fitness = W · similarity to mother − length (default W = 0, i.e. just shortness) |

## The problem

In typogenetics, a **strand** is a string over `A C G T`. A strand is read two bases at a time, and each pair is translated into a *typo-acid* (`cut`, `del`, `swi`, `mvr`, `cop`, `ina`, …), with `AA` acting as a gene separator. So a strand encodes a list of **typo-enzymes**, little programs that run on a strand: they bind to it, move along it, copy it, cut it, insert and delete bases. The strand is the source code and the enzyme is the machine code, which can be fed source code as input.

Hofstadter's challenge (GEB, pp. 512–513): find a strand that, when its own enzymes act on it, and the resulting "daughter strands" are translated and act on themselves, and so on, eventually yields **two copies of the original strand**. A typogenetic quine, more or less.

<!-- Since a quine exists in every Turing-complete language, a proof that typogenetics is Turing complete would also prove that self-replicators exist. As far as I know no such proof has been published. -->

### Previous work

- **Snare (1999)**, an honours thesis, tried the same kind of automated search as this project and found nothing. Most likely because, for simplicity, he made enzyme binding deterministic (see below).
- **Morris (1989)** gave a method for *constructing* self-replicators by hand. Snare adapted it to his interpretation of the rules, and it applies to ours too, so we knew in advance that what we were looking for exists.


## Where this implementation pins down GEB's rules

GEB's description is informal, so some choices had to be made. Most of the rules are implemented as described in the book (and documented in `typogenetics.py`); these are the points worth knowing about.

**Binding is nondeterministic.** An enzyme binds to *any* occurrence of its binding preference in the strand. Every choice is explored, so self-applying a strand produces several possible outcomes, which I call *realities*. This is the main difference from Snare's search, and it is what made self-replicators findable: more possible outcomes means more paths back to the starting strand. (The code also contains deterministic "bind to the rightmost occurrence" variants of `apply` and `daughters`, but the search doesn't use them.) An enzyme consisting of a single typo-acid has no folds, so by convention its binding preference is `A`.

**Molecules and the typo-plasm.** While an enzyme works, its substrate is a double-stranded *typo-molecule*: the primary sequence plus a complementary sequence, either of which may contain gaps. A molecule must stay in one piece. Whenever an operation disconnects it (e.g. `cut`, `del`, or an insertion with copy mode off that leaves a gap opposite a gap), the part the enzyme isn't on drops into the **typo-plasm** and is out of reach for the rest of that enzyme's run. When the enzyme finishes, the substrate and everything in the typo-plasm is broken apart into its single strands (complementary strands read right to left), and that multiset is the enzyme's output.

**Searching in copy mode.** `rpy`/`rpu`/`lpy`/`lpu` may move across gaps in the primary sequence. If copy mode is on, they fill those gaps with complements of the bases above them. Such filled-in bases don't count as search hits.

**Self-application** (following Snare). To self-apply a strand S:

1. Translate S into enzymes e₁, e₂, …, eₙ.
2. Feed S to e₁.
3. Feed each subsequent enzyme the strand its predecessor was bound to when it stopped (the *active strand*).
4. The output is everything all enzymes produced, minus the intermediate active strands consumed as input.

The process stops early if the active strand no longer contains the next enzyme's binding preference, or if it was deleted entirely. Each multiset of outputs (one per reality) is a set of **daughters** of S. Daughters of daughters, and so on, are **descendants**.

**Self-replicator:** S is a self-replicator if, within one consistent reality, its descendants contain at least two copies of S. A **strict** self-replicator is one where both copies appear in the same generation.

## The solution

### Self-replication as a graph problem

Consider the *daughter graph*: a directed multigraph with an edge s → d labelled *i* for each copy of d produced by self-applying s in reality *i*. In it, self-replicators are exactly the strands that lie on **two distinct cycles which agree on realities**. That is, wherever the two cycles share a node, they leave it along edges with the same label (the same outcome of the same nondeterministic choice). Strict self-replicators lie on two such cycles **of equal length**.

The graph is infinite, so it is built incrementally. Each time a strand is *expanded* (self-applied, adding edges to its daughters), every new edge u → v is checked for closing cycles by backtracking over all simple paths from v back to u. Each new cycle is compared with the stored cycles it shares a node with.

### Compressed graph

To keep cycle enumeration cheap, the search builds a *compressed* graph with at most one edge per (source, target) pair. Each edge stores only:

- the set of reality labels it represents, and
- a "duplicate edge" flag that's true iff some reality produced **two or more** copies of the target.

That's enough to decide self-replication:

- Two compressed cycles C₁ and C₂ encode two compatible real cycles iff for every shared node, the label sets on its outgoing edges in C₁ and C₂ intersect. All nodes in C₁ ∩ C₂ are then self-replicators.
- A single compressed cycle containing a duplicate edge already encodes two compatible cycles (take two copies of the duplicated edge, and the same edge everywhere else), so all of its nodes are self-replicators.

On 100,000 random strands of length 0–21, this reduced the number of edges by a factor of about 0.27.

### Search order

**Lexicographic.** Expand strands in order of length, then alphabetically. Simple, and surprisingly effective (see results). Statistics from this search showed that almost all useful cycles are short, while the backtracking spends most of its time deep in the recursion. Limiting cycle length with `-c` therefore saves a lot of time.

**Heuristic.** Put each newly discovered daughter in a priority queue by *fitness*, so the search expands strands that already have many ancestors in the graph and is more likely to close cycles.

- **Shortness** (fitness = −length) alone gets stuck. For example, `GA` has a single daughter `GAA`, whose single daughter is `GAAA`, and so on forever.
- Fixed by **alternating**: expand one fresh strand (the next lexicographic strand not yet seen), then *P* strands from the queue.
- A second fitness term, **similarity to the mother strand**, measured as |LCS(s, mother)| / max(|s|, |mother|), approximates "how few typogenetic operations separate the two". The LCS length is memoized across all strand pairs, using an injective 3-bits-per-base hash of each strand as the cache key. Similarity alone has the same runaway-chain problem, so it is combined with shortness: fitness = −|s| + W · sim(s, mother).

## Results

| search | settings | time limit | expansions | self-replicators | of which strict |
|---|---|---|---|---|---|
| lexicographic | no cycle limit | 1.5 min | 63,000 | 7 | 0 proven |
| lexicographic | `-c 13` | 1.5 min (same 7 found ~4× sooner) | 122,000 | 7 | 2 |
| heuristic, shortness | `-c 20 hrs -p 3` | 30 s | 500 | **38** | **18** |
| heuristic, shortness + similarity | `-c 20 hrs -p 2 -s 15` | 1 min | 880 | **49** | **41** |

The heuristic runs are time-limited, so exact counts vary slightly with machine speed.

The first self-replicators found by lexicographic search were `GAGGAGG`, `GAGGAGGA`, `AGGGGA`, `GGAA`, `GGAAG`, `GGAAGG` and `GGAGGAGG`, all lying on these two cycles (edge labels are reality numbers):

```
GAGGAGG -4-> GAGGAGGA -2-> AAGGA -1-> AAGGGA -0-> AGAAGGGA -3-> AGAGGGA -0->
  AGGGGA -2-> GGAA -1-> GGAAG -1-> GGAAGG -0-> GGAGGAGG -0-> (back to GAGGAGG)

GAGGAGG -4-> GAGGAGGA -2-> GAGGA -0-> GAAGGGA -0->
  AGGGGA -2-> GGAA -1-> GGAAG -1-> GGAAGG -0-> GGAGGAGG -0-> (back to GAGGAGG)
```

The cycles have lengths 11 and 9, so they don't prove strictness. With a cycle limit of 13, `GAGGAGG` and `GAGGAGGA` were later proven strict via a second cycle of length 9.

In the heuristic searches, a further five minutes after the first produced very little new. Run with `-vv` to print the cycles proving each find.

### Why only A and G?

**Every self-replicator found consists only of the bases A and G.** The reason is in the typogenetic code. The four codons made of A and G translate to the gene separator (`AA`), `del` (`AG`), `ina` (`GA`) and `ing` (`GG`). Self-applying an A/G strand therefore only ever produces shorter or longer A/G strands, which makes that region of the daughter graph very densely connected and full of cycles. Hofstadter quite likely designed the code this way on purpose, as a hint for attentive readers hunting for self-replicators by hand.

## Ideas for further work

- **Longer cycles.** Raise the cycle-length limit by maintaining an index of paths between high-out-degree nodes, as in Alibaba's GraphS system for real-time constrained cycle detection (Qiu et al., 2018).
- **Better fitness functions.**
  - Use edit distance with typo-acid-like operations<!--, possibly via the Masek–Paterson algorithm. -->
  - Count self-application steps. Weighting edges this way would allow searching for the *fastest* replicators (cheapest cycles).
  - Inherit fitness from ancestors further back.
  - Reward strands in dense regions of the graph.
  - Use the enzymes a strand translates to.
- **More nondeterminism** (enzyme execution order, choice of input strand) might make self-replicators common enough to *learn* a predictor of self-replication and use it as fitness.
- **Probability of self-replication.** Given a distribution over nondeterministic outcomes, edges get probabilities, and one could compute (or lower-bound, from a subset of descendants) the probability that a strand replicates.

## References

- D. R. Hofstadter, *Gödel, Escher, Bach: an Eternal Golden Braid*, 20th anniversary ed., Basic Books, 1999.
- H. C. Morris, "Typogenetics: A logic for artificial life", *Artificial Life*, 1989, pp. 369–395.
<!-- - L. Varetto, "Typogenetics: An artificial genetic system", *Journal of Theoretical Biology* 160(2), 1993. -->
- A. Snare, *Typogenetics*, honours thesis, Monash University, 1999.
<!-- - V. Kvasnicka, J. Pospichal, "Autoreplicators and hypercycles in typogenetics", *Journal of Molecular Structure: THEOCHEM* 547, 2001. -->
- X. Qiu et al., "Real-time constrained cycle detection in large dynamic graphs", *PVLDB* 11(12), 2018.
<!-- - W. J. Masek, M. S. Paterson, "A faster algorithm computing string edit distances", *JCSS* 20(1), 1980. -->