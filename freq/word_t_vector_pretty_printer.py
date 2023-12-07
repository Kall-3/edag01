import gdb

class PrintWordVector(gdb.Command):
    """Prints a vector of word_t structs until a null string is encountered."""

    def __init__(self):
        print("PrintWordVector loaded")
        super(PrintWordVector, self).__init__("print_word_vector", gdb.COMMAND_USER)

    def invoke(self, arg, from_tty):
        # Split the argument to get the address of the vector and its size.
        args = gdb.string_to_argv(arg)
        if len(args) != 1:
            raise gdb.GdbError("Usage: print_word_vector <address>")

        try:
            address = gdb.parse_and_eval(args[0])
            idx = 0;
    
            # Read and print word_t structs until a null string is encountered.
            while True:
                try:
                    word_ptr = address['word']
                    word = word_ptr.string()
                except gdb.error:
                    raise gdb.GdbError("Invalid memory address for 'word'.\n")

    
                # Stop if the word is a null string.
                if word == '':
                    break
    
                freq = address['freq']
                print(f"{idx}: Word: {word}, Frequency: {freq}")
    
                # Move to the next struct in the vector.
                address = address + 1
                idx = idx + 1
        except gdb.error as e:
            print(f"Error: {e}\n")
    
PrintWordVector()
