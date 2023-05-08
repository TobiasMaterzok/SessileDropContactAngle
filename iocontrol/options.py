import sys

def get_options(argvs, chars, datatypes, defaults):
    """
    Parse command line arguments and return corresponding values based on specified datatypes and default values.

    Args:
    - argvs (list): list of command line arguments.
    - chars (list): list of single character strings that define the valid options for the command line arguments.
    - datatypes (list): list of strings that define the expected data types for each option.
    - defaults (list): list of default values for each option.

    Returns:
    - values (list): list of values corresponding to the command line arguments.
    """
    values = defaults
    for argv in argvs[1::2]:
        if argv in chars:
            datatype = datatypes[chars.index(argv)]
            if datatype == 'int':
                values[chars.index(argv)] = int(argvs[argvs.index(argv) + 1])
            elif datatype == 'float':
                values[chars.index(argv)] = float(argvs[argvs.index(argv) + 1])
            elif datatype == 'str':
                values[chars.index(argv)] = argvs[argvs.index(argv) + 1]
            else:
                print('A data type given in the datatype array, ' + datatype + ' is not valid.')
                sys.exit()
        else:
            print('The option given in the command line, ' + argv + ' is not valid.')
            options = '['
            for op in chars:
                options += op
                if op != chars[-1]:
                    options += ', '
            options += ']'
            print('The possible options are ' + options + '.')
            sys.exit()
    return values
