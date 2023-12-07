define disp
    source word_t_pretty_printer.py
    source word_t_vector_pretty_printer.py
    set print addr off
    tui layout src
    break main
    run < input > output
    display found
    display line_nbr
    display line
end

define n
    next
    print_word_vector words
end
