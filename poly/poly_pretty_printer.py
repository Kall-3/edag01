import gdb

class PrintPolynom(gdb.Command):
    """Prints a polynomial represented by a poly_t vector."""

    def __init__(self):
        super(PrintPolynom, self).__init__("print_polynom", gdb.COMMAND_USER)
        print("PrintPolynom loaded")

    def format_term(self, coef, exp):
        term = ""

        # Handle the coefficient
        if coef < 0:
            term += "- "
        elif i > 0:
            term += "+ "

        coef = abs(coef)
        
        if abs(coef) != 1 or exp == 0:
            term += f"{coef}"

        # Handle 'x' and exponent
        if exp != 0:
            term += "x"
            if exp != 1:
                term += f"^{exp}"

        return term

    def invoke(self, arg, from_tty):
        args = gdb.string_to_argv(arg)
        if len(args) != 1:
            raise gdb.GdbError("Usage: print_polynom <address>")

        try:
            address = gdb.parse_and_eval(args[0])

            terms = []
            i = 0
            for _ in range(10000004):
                element = address.cast(gdb.lookup_type('poly_t').pointer())[i]
                coef = int(element['coef'])
                exp = int(element['exp'])
                if (coef != 0 and exp != 0):
                    formatted_term = self.format_term(coef, exp)
                    terms.append(formatted_term)

            print(" ".join(terms))
        
        except gdb.error as e:
            print(f"Error: {e}")

PrintPolynom()

