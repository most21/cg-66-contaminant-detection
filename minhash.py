import heapq
import sys
import xxhash

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

class MinHash:
    """
        Class to represent MinHash sketches.
        This code is adpated from my solution to HW3 in EN.601.446/646.Sp21
    """

    def __init__(self, data, sketch_size, seed=0):
        """
            Construct a MinHash sketch.
            Args:
                data - list of strings (tokens) to be hashed
                sketch_size - how many hash values to maintain in each sketch
                seed - optionally, set the random seed used in the hash function
        """
        # Hash function
        self._hasher = HashXX32(seed)
        self._max_hash_val = 2**32

        # Multiply every item in the heapq data structure by -1 to form a max heap.
        # Thereby, self._heapq[0] * -1 = max([i * -1 for i in self._heapq])
        self._heap = [-self._max_hash_val for i in range(sketch_size)]
        heapq.heapify(self._heap)

        # Maintain a set of tokens we're currently storing so we do not hash the same token twice
        self._current_tokens = set()

        # Create MinHash sketch from the provided data
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

if __name__ == "__main__":
    sys.path.insert(1, '/Users/matthewost/Documents/Graduate School/(1) Fall 2021/Computational Genomics/Final Project/cg-66-contaminant-detection/data_utils')

    f = "./test_data/fasta_test.fa"
    from data_utils import read_fasta_files, read_fastq_files

    FA = read_fasta_files("./test_data/phix.fa")
    FQ = read_fastq_files("./test_data/phix_reads20.fastq")
    # print(FQ)
    # for r in FQ:
    #     print(r)
    # quit()

    k = 3
    fasta_kmers = [FA[i:i + k] for i in range(len(FA) - k + 1)]
    #print(fasta_kmers)
    print(len(fasta_kmers))
    #quit()
    ref_sketch = MinHash(fasta_kmers, sketch_size=5)

    read_sketches = []
    for read in FQ:
        read_kmers = [read["seq"][i:i + k] for i in range(len(read["seq"]) - k + 1)]
        read_sketches.append(MinHash(read_kmers, sketch_size=5))

    # print(len(read_sketches))
    # test = read_sketches[-1]
    # print(test._heapq)
    # print(test._max_hash_val)
    # quit()

    for i, read_sketch in enumerate(read_sketches):
        print(i, "Sim = ", ref_sketch.similarity(read_sketch), "    Containment = ", read_sketch.containment(ref_sketch))
    pass
