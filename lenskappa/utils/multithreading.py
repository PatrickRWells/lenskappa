import logging
import multiprocessing as mp

class MultiThreadObject:
    num_threads = 1

    @classmethod
    def set_num_threads(cls, num_threads, *args, **kwargs):
        cls.num_threads = num_threads
        cls.validate_num_threads()

    @classmethod
    def get_num_threads(cls, *args, **kwargs):
        return cls.num_threads
    
    @classmethod
    def validate_num_threads(cls, *args, **kwargs):
        num_cores = mp.cpu_count()
        if num_cores < cls.num_threads:
            logging.warning("You requested more cores than this machine has available. "\
                            "I will reduce the number of threads to {} (down from {})".format(num_cores, cls.num_threads))
            cls.set_num_threads(num_cores)

    
    def __init__(self, init_function, *args, **kwargs):
        """
        Wrapper class for handling objects that need to be instantiated seperately in different threads.
        Original use case is np.random.uniform, which is not thread-safe.        
        """
        self._fn = init_function
        self.initialize()

    def initialize(self, *args, **kwargs):
        self.validate_num_threads()
        self._objs = [self._fn() for _ in range(self.num_threads)]

    def __getitem__(self, index):
        if len(self._objs) != self.num_threads:
            self.initialize()
        if index >= self.num_threads:
            logging.error("Unable to get object for thread {}, the total number of threads is only {}".format(index, self.num_threads))
            return
        return self._objs[index]