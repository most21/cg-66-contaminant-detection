from data_utils.dataset import Dataset

class FASTA(Dataset):
    def __init__(self, filename):
        """ Construct a FASTA Dataset object using the data in the given file. """
        # Check that filename has a proper FASTA extension
        assert filename.endswith(".fasta") or filename.endswith(".fa"), filename

        self.filename = filename
        self.data, self.ids = FASTA._read_file(filename)
        self.sizes = tuple([len(frag) for frag in self.data])

    @staticmethod
    def _read_file(filename):
        """
            Read the contents of a given FASTA (or multi-FASTA) file.
            Returns a list of strings with each fragment in the file.
            This function is adapted from one I wrote to read the unitig file in HW4.
        """
        # Read and clean all lines from the file
        with open(filename, "r") as f:
            lines = [l.strip() for l in f.readlines()]

        names = [lines[0][1:]] # initialize with the name of the first fragment
        fragments = []
        cur_frag = ""
        for i in range(1, len(lines)):
            line = lines[i]
            if line[0] == ">": # if we've reached a new fragment
                fragments.append(cur_frag)
                cur_frag = ""
                names.append(line[1:])
            else: # else, just add on the line to the current fragment
                cur_frag += line

        # Append final fragment
        fragments.append(cur_frag.upper())

        return fragments, names

    def __str__(self):
        """ Generate string summary representation of FASTA dataset. """
        output = "FASTA Dataset\n"
        output += "Filename: " + self.filename + "\n"
        output += "IDs: " + str(self.ids) + "\n"
        output += "Sizes: " + str(self.sizes)
        return output

    def __iter__(self):
        """ Custom iterator over FASTA classes that just iterates over the fragments. """
        return iter(self.data)

    def __getitem__(self, idx):
        """ Allows integer indexing of FASTA objects- return the value in a certain position in the data string. """
        return self.data[idx]

    def __len__(self):
        """ Return the number of fragments in the FASTA data. """
        return len(self.sizes)

    def get_data(self):
        """ Return the raw FASTA data (list of strings). """
        return self.data
