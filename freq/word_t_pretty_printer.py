import gdb

class WordTPrettyPrinter:
    def __init__(self, val):
        self.val = val

    def to_string(self):
        # Access struct fields
        word = self.val['word'].string()
        freq = self.val['freq']
        return "word = {}, freq = {}".format(word, freq)

def lookup_type(val):
    # Adjust the type name based on your program's type name (including namespaces if in C++)
    if str(val.type) == "word_t":
        return WordTPrettyPrinter(val)
    return None

gdb.pretty_printers.append(lookup_type)
