import heapq
import numpy as np
import sys
import xxhash

from data_utils.data_utils import save_object_to_file, load_object_from_file

CACHE_PATH = "cache/minhash/"

class HashXX32:
    """
        Class to represent a hash function.
        This Hash class is taken from the HW3 starter code in EN.601.446/646.Sp21
    """
    def __init__(self, seed):
        self.seed = seed
        self.hasher = xxhash.xxh32(seed=seed)

    def hash(self, item):
        self.hasher.reset()
        self.hasher.update(item)
        return self.hasher.intdigest() % sys.maxsize

class HashXX64:
    """
        Class to represent a hash function.
        This Hash class is taken from the HW3 starter code in EN.601.446/646.Sp21
    """
    def __init__(self, seed):
        self.seed = seed
        self.hasher = xxhash.xxh64(seed=seed)

    def hash(self, item):
        self.hasher.reset()
        self.hasher.update(item)
        return self.hasher.intdigest() % sys.maxsize

class MinHash:
    """
        Class to represent MinHash sketches.
        This code is adpated from my solution to HW3 in EN.601.446/646.Sp21
    """

    def __init__(self, sketch_size, data=None, seed=0, id=None):
        """
            Construct a MinHash sketch.
            Args:
                data - list of strings (tokens) to be hashed
                sketch_size - how many hash values to maintain in each sketch
                seed - optionally, set the random seed used in the hash function
        """
        # Id attribute for easier book-keeping
        self.id = id

        # Hash function
        self._hasher = HashXX64(seed)
        self._max_hash_val = 2**64

        # Multiply every item in the heapq data structure by -1 to form a max heap.
        # Thereby, self._heap[0] * -1 = max([i * -1 for i in self._heap])
        self._heap = [-self._max_hash_val for i in range(sketch_size)]
        heapq.heapify(self._heap)

        # Maintain a set of tokens we're currently storing so we do not hash the same token twice
        self._current_tokens = set()

        # Create MinHash sketch from the provided data, if applicable
        if data:
            for token in data:
                self.update(token)

    def update(self, obj):
        """
        Add obj to MinHash sketch.

        We store -1 * (hashed result) in `self._heap`, by which
        (-1 * self._heap[0]) is the maximal value of the sketch.
        """
        h = self._hasher.hash(obj)
        if h < -self._heap[0] and obj not in self._current_tokens:
            # Add new object to the heap
            ret = heapq.heappushpop(self._heap, -h)
            heapq.heapify(self._heap)

            # Keep track of this new object in our set so we can check for duplicates later (i.e. mark it as seen)
            self._current_tokens.add(obj)

            # Remove the value that got removed upon this update from our duplicates set (i.e. mark it as unseen)
            for val in self._current_tokens:
                if self._hasher.hash(val) == -ret:
                    self._current_tokens.remove(val)
                    break

    def get_heap(self):
        """ Return the MinHash sketch. """
        return self._heap

    def cardinality(self):
        """ Return the estimated cardinality of the sketch. """
        m_k = -self._heap[0]
        return (len(self._heap) * self._max_hash_val / m_k) - 1

    def union(self, b):
        """ Return the union sketch (A and B). """
        size = len(self.get_heap())

        # Merge two heaps and create new heap
        merged = heapq.merge(self.get_heap(), b.get_heap())
        full_union = list(set(merged)) # remove duplicates in the union
        heapq.heapify(full_union)

        # Extract the n values we care about
        union_sketch = heapq.nlargest(size, full_union)
        heapq.heapify(union_sketch)

        return union_sketch

    def size_union(self, B):
        """ Return the size of union sketch (|A and B|). """
        union_sketch = self.union(B)
        m_k = -union_sketch[0]
        return (len(self._heap) * self._max_hash_val / m_k) - 1

    def containment(self, B):
        """
            Estimate the containment of A in B.
            Calculate: |sketch(A).intersection(sketch(B))| / |sketch(A)|
        """
        size_A = self.cardinality()
        size_B = B.cardinality()
        size_AB_union = self.size_union(B)
        return (size_A + size_B - size_AB_union) / size_A

    def similarity(self, B):
        """
            Estimate the Jaccard similarity between A and B.
            Calculate: |sketch(A).intersection(sketch(B))| / |sketch(A).union(sketch(B))|
        """
        size_A = self.cardinality()
        size_B = B.cardinality()
        size_AB_union = self.size_union(B)
        return (size_A + size_B - size_AB_union) / size_AB_union


def make_kmers(dna, k):
    """ Generate list of k-mers from the given DNA string.  """
    return [dna[i:i + k] for i in range(len(dna) - k + 1)]

# def containment(A, B):
#     size_A = A.count()
#     size_B = B.count()
#
#     union = A.copy()
#     union.merge(B)
#     size_union = union.count()
#
#     #print(size_A, size_B, size_union)
#     return (size_A + size_B - size_union) / size_A

def build_ref_sketches(fasta_objects, k, sketch_size, cache=True):
    """
        Given a list of FASTA objects, create a MinHash sketch for each one (or load existing sketch from cache).
        Args:
            k - int, the k-mer size
            sketch_size - int, the number of hashes to store in each MinHash sketch
    """
    ref_sketch_ids = []
    ref_sketches = []
    for fasta in fasta_objects:
        sketch_id = "//".join(fasta.ids)
        filename = fasta.filename.split("/")[-1]

        # Try loading a cached ref sketch, if it exists. Else we will construct it from scratch.
        try:
            sketch = load_object_from_file(CACHE_PATH + filename + ".b")
            ref_sketches.append(sketch)
            ref_sketch_ids.append(sketch_id)
            continue
        except FileNotFoundError:
            print(f"[WARNING] Could not find cached sketch for '{filename}'. Building from scratch...")

        # Break each fragment in the FASTA file into chunks and each chunk into k-mers and add it to our MinHash sketch
        #sketch = MinHash(num_perm=256)
        sketch = MinHash(sketch_size)
        for frag in fasta:
            # Break each fragment into kmers and create sketch
            # For memory efficiency, we generate each k-mer one at a time and update the sketch in a stream
            for i in range(len(frag) - k + 1):
                kmer = frag[i: i + k]
                sketch.update(kmer)

        ref_sketches.append(sketch)
        ref_sketch_ids.append(sketch_id)

        # If cache is set to True, save this sketch to disk
        # Note that the _hasher attribute must be set to None to play nice with pickle.
        if cache:
            sketch._hasher = None
            save_object_to_file(sketch, CACHE_PATH + filename + ".b")

    return ref_sketches, ref_sketch_ids


def classify_reads(fastq_obj, cont_ref_sketches, cont_ref_sketch_ids, des_ref_sketches, des_ref_sketch_ids, k, sketch_size):
    """ Label each read as contaminated, desired, or other. """
    results = {"Contaminated": [], "Desired": [], "Unassigned": []}

    # Create MinHash sketch for each read in the file using its k-mers
    for i, read in enumerate(fastq_obj):
        # Create empty sketch for this read
        read_sketch = MinHash(sketch_size)

        # Break read into k-mers and populate the sketch
        kmers = make_kmers(read["seq"], k)
        for kmer in kmers:
            read_sketch.update(kmer)

        # Compare read sketch to each reference sketch and compute mean score
        mean_cont_score = np.max([read_sketch.containment(cont_ref_sketch) for cont_ref_sketch in cont_ref_sketches])
        mean_des_score = np.max([read_sketch.containment(des_ref_sketch) for des_ref_sketch in des_ref_sketches])

        # Compute max score (which reference this read likely came from)
        if mean_cont_score > mean_des_score:
            results["Contaminated"].append(read)
        elif mean_cont_score < mean_des_score:
            results["Desired"].append(read)
        else:
            results["Unassigned"].append(read)

    return results


def minhash_engine(cont_fasta_objs, des_fasta_objs, fastq_obj, k=21, sketch_size=1000):
    # Build sketches for each reference genome, contaminants and desired
    cont_ref_sketches, cont_ref_sketch_ids = build_ref_sketches(cont_fasta_objs, k, sketch_size, cache=True)
    des_ref_sketches, des_ref_sketch_ids = build_ref_sketches(des_fasta_objs, k, sketch_size, cache=True)

    # For each sequencing read in each FASTQ file, identify whether it is a contaminant or not
    results = classify_reads(
        fastq_obj,
        cont_ref_sketches,
        cont_ref_sketch_ids,
        des_ref_sketches,
        des_ref_sketch_ids,
        k,
        sketch_size)

    return results
