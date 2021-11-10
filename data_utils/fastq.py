from dataset import Dataset

class FASTQ(Dataset):

    ######### Nested class to define an iterator over FASTQ objects ############
    class FASTQIterator:
        """ Custom iterator for FASTQ files that can iterate over any of the 3 major series (or all of them). """
        def __init__(self, dataset, series="all"):
            series_dict = {
                "id": dataset.read_ids,
                "seq": dataset.read_data,
                "quality": dataset.read_qualities,
                "all": list(zip(dataset.read_ids, dataset.read_data, dataset.read_qualities))
                }
            self._dataset = dataset
            self._idx = 0
            self._series_label = series
            self._series = series_dict[series]

        def __next__(self):
            """ Return next value in the specified series. """
            if self._idx < len(self._series):
                result = self._series[self._idx]
                self._idx += 1

                # If iterating over all 3 series, construct a dictionary return value
                if self._series_label == "all":
                    n, r, q = result
                    result = {"id": n, "seq": r, "quality": q}
                return result

            # Once iteration is finished, reset iterator series to be "all" by default
            self._dataset._iterator_series = "all"
            raise StopIteration
    ################## End of FASTQIterator nested class ###################

    def __init__(self, filename, id=None):
        """ Construct a FASTQ Dataset object using the data in the given file. """
        # Check that filename has a FASTQ extension
        assert filename.endswith(".fastq")

        self.filename = filename
        self.id = id if id else filename
        self.read_ids, self.read_data, self.read_qualities = FASTQ._read_file(filename)
        self.num_reads = len(self.read_ids)
        self._iterator_series = "all"

    @staticmethod
    def _read_file(filename):
        """
            Read the contents of a given FASTQ file. Returns 3 lists- IDs, DNA seqs, and quality vals.
            This function is adapted from the provided code for HW2.
        """
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
            # Find the read with this id - O(n)
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

    def __len__(self):
        """ Return the number of reads in this FASTQ dataset. """
        return self.num_reads

    def get_read_ids(self):
        """ Return list of read IDs. """
        return self.read_ids

    def get_read_sequences(self):
        """ Return list of read sequences. """
        return self.read_data

    def get_read_qualities(self):
        """ Return list of read quality values. """
        return self.read_qualities
