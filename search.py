"""search for self-replicators"""

from __future__ import annotations
from typogenetics import *
import itertools
import heapq
from collections import defaultdict
import pickle
from datetime import datetime
import os
from matplotlib import pyplot as plt
from functools import cached_property
import numpy as np
import argparse


def random_strand(length: int = None) -> Strand:
    if length is None:
        length = np.random.randint(Strand.MAX_64BIT_LEN + 1)

    return Strand([Base(np.random.randint(len(Base))) for _ in range(length)])


def generate_strands(min_len=0, max_len: int = None) -> Iterator[Strand]:
    """generate all strands in given length range, in lexicographic order"""
    length = min_len

    while max_len is None or length <= max_len:
        for seq in itertools.product(Base, repeat=length):
            yield Strand(bases=list(seq))

        length += 1



@dataclass
class Daughter:
    """edge in the compressed family tree"""

    strand: Strand
    """the daughter strand"""

    reality2num_copies: Dict[int,int]
    """mapping from reality ID to the number of copies of this
    strand \"conceived\" in this reality"""

    @cached_property
    def realities(self) -> Set[int]:
        """the IDs of the realities in which this daughter was conceived"""
        return set(self.reality2num_copies.keys())

    @cached_property
    def duplicate_reality(self) -> int | None:
        """one of the realities in which several copies were conceived (if any)"""
        try:
            return next(rlty for rlty, n in self.reality2num_copies.items() if n >= 2)
        except StopIteration:
            return None

    def __repr__(self) -> str:
        rlts = (str(r) + (f"x{n}" if n > 1 else "")
                for r, n in self.reality2num_copies.items())

        return f"<{self.strand} {','.join(rlts)}>"



def compressed_daughters(mother_strand: Strand) -> Iterator[Daughter]:
    
    strand2reality_dict: Dict[Strand, Dict[int,int]] \
        = defaultdict(lambda: defaultdict(int))
        # strand -> (realityID -> #copies in that reality)

    for realityID, daughters in enumerate(mother_strand.daughter_sets):
        for dtr_strand in daughters:
            strand2reality_dict[dtr_strand][realityID] += 1

    for strand, rlty2n_copies in strand2reality_dict.items():
        yield Daughter(strand, rlty2n_copies)



Cycle = List[Daughter]
"""Cycle in the compressed family tree.
Thee first daughter's edge goes from end to start of list"""

def realities_from(c: Cycle) -> dict[Strand, Set[int]]:
    """mapping strand -> realities in which the cycle successor is its daughter"""
    return {c[i - 1].strand: c[i].realities for i in range(len(c))}


def contains_several_cycles(compressed_cycle: Cycle) -> bool:
    """does it encode several cycles in the non-compressed family tree?"""
    stats["cycles of length"][len(compressed_cycle)] += 1

    return any(edge.duplicate_reality is not None for edge in compressed_cycle)



def new_cycles(daughters_of: Dict[Strand, List[Daughter]], source: Strand,
               new_edge: Daughter, max_cycle_len: int = None) -> Iterator[Cycle]:

    """yields the simple cycles incurred by the new edge"""

    visited: Set[Strand] = set() 
    path: List[Daughter] = []

    def simple_paths_from(edge: Daughter, depth=0) -> Iterator[Cycle]:
        """paths from this edge's destination to the source of the new edge"""

        stats["calls at depth"][depth] += 1

        visited.add(edge.strand)
        path.append(edge)

        if edge.strand == source:
            yield path.copy()
        
        elif max_cycle_len is None or depth < max_cycle_len - 1:
            for daughter in daughters_of[edge.strand]:
                if daughter.strand not in visited:
                    yield from simple_paths_from(daughter, depth + 1)
        
        path.pop()
        visited.remove(edge.strand)

    yield from simple_paths_from(new_edge)


def self_reps_on(c1: Cycle, c2: Cycle) -> SelfRepBatch | None:
    """the self-replicators common to these cycles,
    or None if the cycles don't agree on reality"""

    c1_realities_from = realities_from(c1)
    c2_realities_from = realities_from(c2)

    for strand in c1_realities_from:
        if strand in c2_realities_from and \
            not c1_realities_from[strand].intersection(c2_realities_from[strand]):
                return None

    return SelfRepBatch(c1, c2)



class StrandQueue:
    """min-heap of strands"""

    heap: List[Strand]
    strands_inside: Set[Strand]

    def __init__(self):
        self.heap = []
        self.strands_inside = set()

    def push(self, s: Strand, priority: float) -> None:
        heapq.heappush(self.heap, (priority, s))
        self.strands_inside.add(s)

    def pop(self) -> Strand | None:
        try:
            _, s = heapq.heappop(self.heap)
            self.strands_inside.remove(s)
            return s
        except IndexError:
            return None

    # OPT use heappushpop

    def __contains__(self, s: Strand) -> bool:
        return s in self.strands_inside

    def __len__(self) -> int:
        return len(self.heap)



class SearchState:
    daughters_of: Dict[Strand, List[Daughter]] 
    """the compressed daughter relation graph"""

    cycles_containing: Dict[Strand, List[Cycle]]
    """map strand -> the cycles it's a part of"""

    queue: StrandQueue
    """the strands awaiting expansion"""

    def __init__(self, fname: str = None):
        if fname is None:
            self.daughters_of = defaultdict(list)
            self.cycles_containing = defaultdict(list)
            self.queue = StrandQueue()
        else:
            with open(fname, "rb") as f:
                members = pickle.load(f)

            self.__dict__.update(members)


    def save(self) -> None:
        stamp = datetime.now().strftime(r"%d-%m@%H-%M-%S")
        fname = os.path.expanduser("~") + "\\src\\search_state\\" + stamp

        with open(fname + ".pickle", mode="xb") as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)



class SelfRepBatch:
    """a batch of unique self-replicators"""

    strands: List[Strand]
    cycle1: Cycle
    cycle2: Cycle

    def __init__(self, fst_cycle: Cycle, snd_cycle: Cycle = None):
        self.cycle1 = fst_cycle

        if snd_cycle is None:
            self.strands = [edge.strand for edge in fst_cycle]
            self.cycle2 = None
        else:
            fst_strands = set(edge.strand for edge in fst_cycle)
            self.cycle2 = snd_cycle
            self.strands = [edge.strand for edge in snd_cycle \
                                         if edge.strand in fst_strands]

        assert len(self.strands) == len(set(self.strands))


    def __repr__(self) -> str:
        return repr(self.strands)

    def __iter__(self) -> Iterator[Strand]:
        return iter(self.strands)

    def __len__(self) -> int:
        return len(self.strands)

    @property
    def cycles(self) -> Iterator[Cycle]:
        yield self.cycle1
        if (c2 := self.cycle2) is not None: yield c2

    @property
    def strict(self) -> bool:
        """do copies of these strands appear in the same generation?"""
        return self.cycle2 is None or len(self.cycle1) == len(self.cycle2)
    
    
    def prove(self) -> None:
        """print the two distinct (un-compressed) cycles"""
        if self.cycle2 is None:
            print(f"[len={len(self.cycle1)}]", end=" ")

            double_edge = False

            for edge in self.cycle1:
                if not double_edge and \
                    (dup := edge.duplicate_reality) is not None:
                    print(f" ={dup}=> ", end="")
                    double_edge = True
                else:
                    print(f" -{next(iter(edge.realities))}-> ", end="")

                print(edge.strand)
            
            print()
        else:
            out_realities = [realities_from(c) for c in self.cycles]
            common_reality: Dict[Strand, int] = {}

            for strand in (c1_OR := out_realities[0]):
                if strand in (c2_OR := out_realities[1]):
                    intersection = c1_OR[strand].intersection(c2_OR[strand])
                    common_reality[strand] = next(iter(intersection))
            
            cycle_len = [len(c) for c in self.cycles]

            for i, cycle in enumerate(self.cycles):
                print(f"[len={cycle_len[i]}]", end=" ")
                for edge in cycle:
                    try:
                        r = common_reality[s := edge.strand]
                    except KeyError:
                        r = next(iter(out_realities[i][s]))

                    print(f"{edge.strand} -{r}-> ", end="")
                
                print()


FitnessFunction = Callable[[Strand,Strand], float]
"""(a strand, its mother) -> the strand's fitness value"""

def expand(mother_strand: Strand, state: SearchState,
           cycle_len_limit: int = None,
           daughter_fitness: FitnessFunction = None) -> Iterator[SelfRepBatch]:
    
    """Yields any found batches of self-reps.
    Enqueues daughters if fitness function is given."""
    
    if args.verbosity >= 3: print("expanding", mother_strand)
    stats["strands of length"][len(mother_strand)] += 1

    for daughter in compressed_daughters(mother_strand):
        state.daughters_of[mother_strand].append(daughter)

        if daughter_fitness is not None and \
            daughter.strand not in state.daughters_of and \
            daughter.strand not in state.queue:
            
            state.queue.push(
                daughter.strand,
                priority=-daughter_fitness(daughter.strand, mother_strand
            ))

        cycle_gen = new_cycles(
            state.daughters_of, source=mother_strand, new_edge=daughter,
            max_cycle_len=cycle_len_limit
        )
        # check if new edge makes cycles
        for new_cycle in cycle_gen:
            if contains_several_cycles(new_cycle):
                yield SelfRepBatch(new_cycle)

            # to avoid checking against the same old cycle multiple times
            compared_with: List[Cycle] = [] 
            
            # check if any strands in this cycle are in other cycles
            for graph_edge in new_cycle:
                for old_cycle in \
                    (old_cycles := state.cycles_containing[graph_edge.strand]):
                    
                    if old_cycle not in compared_with:
                        compared_with.append(old_cycle)
                        if (reps := self_reps_on(old_cycle, new_cycle)) is not None:
                            yield reps

                old_cycles.append(new_cycle) 



def lexicographic_search(cycle_maxlen: int = None) -> Iterator[SelfRepBatch]:
    """yields batches of self-reps (may contain duplicates across batches)"""

    state = SearchState()
    try:
        for strand in generate_strands():
            yield from expand(strand, state, cycle_len_limit=cycle_maxlen)
    except KeyboardInterrupt:
        # state.save()
        return



### fitness functions ###

def shortness(s: Strand, mother: Strand) -> float:
    return -len(s)


LCS_cache: Dict[Tuple[int,int], int] = {}

def LCS_len_cached(h1: int, h2: int) -> int:
    for arg_perm in [(h1, h2), (h2, h1)]:
        try: return LCS_cache[arg_perm]
        except KeyError: pass

    if h1 == 0 or h2 == 0:
        ret = 0
    elif h1 & 0b111 == h2 & 0b111:
        ret = 1 + LCS_len_cached(h1 >> 3, h2 >> 3)
    else:
        ret = max(LCS_len_cached(h1 >> 3, h2), LCS_len_cached(h1, h2 >> 3))

    LCS_cache[h1, h2] = ret
    return ret


def LCS_len(x: Strand, y: Strand) -> int:
    """length of the longest common subsequence"""
    # see CLRS 4th ed., section 14.4
    return LCS_len_cached(hash(x), hash(y))


def similarity(s: Strand, mother: Strand) -> float:
    """between 0 and 1 (1 iff strands equal)"""
    if len(mother) == 0: # daughter also empty strand
        return 1
    else:
        return LCS_len(s, mother) / max(len(s), len(mother))


def combination(similarity_weight: float) -> FitnessFunction:
    """returns fitness function combining shortness and similarity"""
    return lambda s, mother: \
        shortness(s, mother) + similarity_weight * similarity(s, mother)


def plot_distribution(fitness_fun: FitnessFunction, sample_sz=1000) -> None:
    vals = []
    for _ in range(sample_sz):
        s = random_strand()
        for d in compressed_daughters(s):
            vals.append(fitness_fun(d.strand, s))

    plt.hist(vals, bins=50)
    plt.show()

##########################


def heuristic_search(fitness: FitnessFunction, cycle_maxlen: int = None,
                     dequeue_batch=1) -> Iterator[SelfRepBatch]:

    """Yields batches of self-reps (may contain duplicates across batches)."""

    def fresh_strand() -> Strand:
        nxt_strand = next(strand_generator)
        while nxt_strand in state.daughters_of or nxt_strand in state.queue:
            nxt_strand = next(strand_generator)
        
        return nxt_strand


    state = SearchState()
    strand_generator = generate_strands()
        
    while True:
        try:
            yield from expand(fresh_strand(), state,
                        cycle_len_limit=cycle_maxlen, daughter_fitness=fitness)

            for _ in range(dequeue_batch): 
                if (nxt_strand := state.queue.pop()) is None:
                    break

                yield from expand(nxt_strand, state,
                    cycle_len_limit=cycle_maxlen, daughter_fitness=fitness)

            stats["queue size"].append(len(state.queue))

        except KeyboardInterrupt:
            # state.save()
            return



def plot_stats(stats, n_figures=2, fname="stats") -> None:
    
    def autolabel(rects, ax):
        """Attach a text label above each bar displaying its height"""
        for rect in rects:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., height + 0.01,
                    str(int(height)), ha='center', va='bottom')


    _, axs = plt.subplots(1, n_figures)

    for i, ax in enumerate(axs):
        match i:
            case 0:
                ax.set_title("calls at depth (path enumeration)")
                ax.bar(*zip(*stats["calls at depth"].items()))
            case 1:
                ax.set_title("#detected cycles by length")
                rects = ax.bar(*zip(*stats["cycles of length"].items()))
                autolabel(rects, ax)
            case 2:
                ax.set_title("#expanded strands by length")
                rects = ax.bar(*zip(*stats["strands of length"].items()))
                autolabel(rects, ax)
            case 3:
                ax.set_title("queue size over expansion batches")
                ax.plot(stats["queue size"])

    plt.show()
    try:
        plt.savefig(fname + ".pdf")
    except:
        print("stats plot not saved")



def print_time_since(t_start: datetime):
    mins, secs = divmod((datetime.now() - t_start).seconds, 60)
    total_expanded = sum(v for v in stats["strands of length"].values())
    print(f"[t={mins}:{secs:02}, expanded={total_expanded}]")



parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-c", "--max_cycle", metavar="L", type=int, help="cycle length limit")
parser.add_argument("-v", "--verbosity", action="count", default=0,
                    help="increase output verbosity (up to 3 times, eg. -vvv)")
parser.add_argument("-p", "--plot_stats", action="store_true",
                    help="show and save plot of search statistics")

subparsers = parser.add_subparsers(help="search methods", required=True, 
                                   dest="method_name")
subparsers.add_parser("lex", help="lexicographic search")

heur = subparsers.add_parser("hrs", help="heuristic search")
heur.add_argument("-p", "--deq_batch", type=int, metavar="P", default=1,
                  help="dequeue P strands for every new one (default: P=1)")

heur.add_argument("-s", "--similarity_weight", type=float, metavar="W", default=0,
                  help="fitness = W * similarity - length (default: W=0)")

args = parser.parse_args()

if args.method_name == "lex":
    batch_gen = lexicographic_search(cycle_maxlen=args.max_cycle)
else:
    fit_fun = shortness if (sw := args.similarity_weight) == 0 \
              else combination(similarity_weight=sw)
    
    batch_gen = heuristic_search(fitness=fit_fun,
        cycle_maxlen=args.max_cycle, dequeue_batch=args.deq_batch)


stats = {
    "strands of length": defaultdict(int),
    "calls at depth": defaultdict(int),
    "cycles of length": defaultdict(int),
    "queue size": []
}
unique_selfreps: Set[Strand] = set()
strict_selfreps: Set[Strand] = set()
longest_useful_cycle = -1

print("starting search, press ctrl+C to stop")
t0 = datetime.now()

for batch in batch_gen:
    new_strands = [s for s in batch if s not in unique_selfreps]
    new_strict = [s for s in batch if s not in strict_selfreps] if batch.strict else []

    if new_strands or new_strict:
        longest_useful_cycle = max(longest_useful_cycle, *map(len, batch.cycles))

    if new_strands or new_strict or batch.cycle2 is None:
        if args.verbosity >= 1:
            print()
            print_time_since(t0)
            
            if batch.strict: print("STRICT", end=" ")
            print(f"self-rep batch (size={len(batch)}), new:", new_strands)

        unique_selfreps.update(new_strands)

        if new_strict:
            strict_selfreps.update(new_strict)
            if args.verbosity >= 1: print("newly proven strict:", new_strict)

        if args.verbosity >= 2:
            print("on cycles:")
            batch.prove()


if unique_selfreps:
    print(f"\nfound {len(unique_selfreps)} unique self-reps: {unique_selfreps}")
    if strict_selfreps:
        print(f"... {len(strict_selfreps)} of which proven strict:", strict_selfreps)
else:
    print("found no self-reps")

print("\nelapsed ", end=""); print_time_since(t0)
if longest_useful_cycle != -1:
    print(f"longest cycle that proved something new: {longest_useful_cycle} edges")

if args.plot_stats:
    plot_stats(stats, n_figures=3 if args.method_name == "lex" else 4)