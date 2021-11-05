from abc import ABC, abstractmethod # abstract base class (abc)

class Dataset(ABC):

    @staticmethod
    @abstractmethod
    def _read_file(filename):
        """ Open and read a file of a particular format (FASTQ/FASTA). Must be overridden by each subclass. """
        raise NotImplementedError(f"_read_file is not implemented for object of type {type(self)}")

    @abstractmethod
    def __iter__(self):
        """ Overloaded to allow iteration through a dataset. """
        raise NotImplementedError(f"__iter__ is not implemented for object of type {type(self)}")

    @abstractmethod
    def __getitem__(self, idx):
        """ Overloaded to support random access via integer indexing. """
        raise NotImplementedError(f"__get__ is not implemented for object of type {type(self)}")

class FASTA(Dataset):
    def __init__(self, filename):
        """ Construct a FASTA Dataset object using the data in the given file. """
        self.filename = filename
        self.data, self.id = self._read_file(filename)
        self.size = len(self.data)

    @staticmethod
    def _read_file(filename):
        """ Read the contents of a given FASTA file. Returns tuple (data, id) where id is the first row of the file. """
        with open(filename, "r") as f:
            lines = f.readlines()
        return "".join([l.strip() for l in lines[1:]]), lines[0][1:].strip()

    def __str__(self):
        """ Generate string summary representation of FASTA dataset. """
        return "FASTA Dataset\n" + "ID: " + self.id + "\nFilename: " + self.filename + "\n" + "Size: " + str(self.size)

    def __iter__(self):
        """ Custom iterator over FASTA classes that just iterates over the data string. """
        return iter(self.data)

    def __getitem__(self, idx):
        """ Allows integer indexing of FASTA objects- returns a certain position in the data string. """
        return self.data[idx]

class FASTQ(Dataset):

    class FASTQIterator:
        """ Custom iterator for FASTQ files that can iterate over any of the 3 major series: names, reads, quality values"""
        def __init__(self, dataset, series="all"):
            series_dict = {
                "id": dataset.read_ids,
                "seq": dataset.read_data,
                "quality": dataset.read_qualities,
                "all": [{"id": n, "seq": r, "quality": q} for n, r, q in zip(dataset.read_ids, dataset.read_data, dataset.read_qualities)]
                }
            self._dataset = dataset
            self._idx = 0
            self._series = series_dict[series]

        def __next__(self):
            """ Return next value in the specified series. """
            if self._idx < len(self._series):
                result = self._series[self._idx]
                self._idx += 1
                return result
            self._dataset._iterator_series = "all" # reset iterator series to be "all" by default
            raise StopIteration

    def __init__(self, filename, id=None):
        """ Construct a FASTQ Dataset object using the data in the given file. """
        self.filename = filename
        self.id = id if id else filename
        self.read_ids, self.read_data, self.read_qualities = self._read_file(filename)
        self.num_reads = len(self.read_ids)
        self._iterator_series = "all"

    @staticmethod
    def _read_file(filename):
        """ Read the contents of a given FASTQ file. Returns 3 lists- IDs, DNA seqs, and quality vals"""
        names, seqs, quals = [], [], []
        with open(filename, "r") as f:
            while True:
                # If next line is blank, we've reached the end of the file
                line = f.readline()
                if len(line) == 0:
                    break

                # Extract the read ID
                name = line[1:].rstrip()

                # Extract the read data
                seq = f.readline().rstrip()

                # Skip a line (+)
                f.readline()

                # Extract the read qualities
                qual = f.readline().rstrip()

                # Update lists
                names.append(name)
                seqs.append(seq)
                quals.append(qual)

        return names, seqs, quals

    def __str__(self):
        """ Generate string summary representation of the FASTQ dataset. """
        return "FASTQ Dataset\n" + "ID: " + self.id + "\nFilename: " + self.filename + "\n" + "Size: " + str(self.num_reads)

    def __getitem__(self, idx):
        """
            Return a read with all accompanying data (ref ID, quality vals) by its integer index or by its ref ID.
            Note that querying by integer index is O(1) but querying by ref ID is O(n) in the worst case.
            If it turns out we need faster lookup by ref ID, I will incorporate a hash table instead.
         """
        if isinstance(idx, int):
            if idx >= self.num_reads:
                raise IndexError("Dataset index out of range")
            return {
                "id": self.read_ids[idx],
                "seq": self.read_data[idx],
                "quality": self.read_qualities[idx]
                }
        elif isinstance(idx, str):
            # Find the read with this id
            for i, name in enumerate(self.read_ids):
                if idx == name:
                    return {
                        "id": self.read_ids[i],
                        "seq": self.read_data[i],
                        "quality": self.read_qualities[i]
                        }
            raise IndexError("Read ID not found in dataset")

        raise TypeError("Index must be type int or str")

    def __iter__(self):
        """ Iterator over the FASTQ data for whichever series is specified (names, reads, qualities). """
        return self.FASTQIterator(self, series=self._iterator_series)

    def __call__(self, s):
        """
            Make FASTQ objects callable to set the iteration mode (i.e. which series can be iterated over).
            Allowable values: id, seq, quality, all
        """
        if s not in ("id", "seq", "quality", "all"):
            raise ValueError("Illegal value for iteration series")
        self._iterator_series = s
        return self

    # def get_read_data_iterator(self):
    #     """ Iterator over the FASTQ read data. """
    #     return FASTQIterator(self, series="seq")
    #
    # def get_read_name_iterator(self):
    #     """ Iterator over the FASTQ read names. """
    #     return FASTQIterator(self, series="id")
    #
    # def get_read_quality_iterator(self):
    #     """ Iterator over the FASTQ read quality values. """
    #     return FASTQIterator(self, series="quality")



if __name__ == "__main__":
    f = "test_data/phix_reads20.fastq"
    FQ = FASTQ(f)

    for i, item in enumerate(FQ):
        print(i, item)
