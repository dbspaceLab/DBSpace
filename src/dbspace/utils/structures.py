from collections import defaultdict


def nestdict():
    """
    Simple, very powerful, structure for deep nested dictionaries
    """

    return defaultdict(nestdict)
