class Suffix:
    def __init__(self, traversed, count, marker, suffix):
        self.traversed = traversed
        self.count = count
        self.marker = marker
        self.suffixes = [suffix]
