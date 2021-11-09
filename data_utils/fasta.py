from dataset import Dataset

class FASTA(Dataset):
    def __init__(self, filename):
        """ Construct a FASTA Dataset object using the data in the given file. """
        # Check that filename has a proper FASTA extension
        assert filename.endswith(".fasta") or filename.endswith(".fa")

        self.filename = filename
        self.data, self.id = FASTA._read_file(filename)
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
        """ Allows integer indexing of FASTA objects- return the value in a certain position in the data string. """
        return self.data[idx]

    def __len__(self):
        """ Return the length of the FASTA data. """
        return self.size

    def get_data(self):
        """ Return raw string of FASTA data. """
        return self.data
