from abc import ABC, abstractmethod # abstract base class (abc)

class Dataset(ABC):
    """ Abstract data type (ADT) for a dataset object. """

    @staticmethod
    @abstractmethod
    def _read_file(filename):
        """ Open and read a file of a particular format. """
        pass

    @abstractmethod
    def __iter__(self):
        """ Create an iterator through the dataset. """
        pass

    @abstractmethod
    def __getitem__(self, idx):
        """ Enable dataset access via indexing. """
        pass

    @abstractmethod
    def __len__(self):
        """ Return size of dataset. """
        pass
