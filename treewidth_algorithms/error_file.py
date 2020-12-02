class Error(Exception):
    """Base class for other exceptions"""
    pass


class MemoryLimitViolation(Error):
    """Raised when too much memory is used"""
    pass

class NoPersistentCycle(Error):
    """Raised when no perisistent cycle is found"""
    pass