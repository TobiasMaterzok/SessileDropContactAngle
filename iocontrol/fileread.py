import numpy as np

def file_to_array(filename, row_len=-1):
    # Define list of comment indicators
    comment_indicators = ['#', '@', '\n']
    # Open file for reading
    file_pointer = open(filename, "r")
    data = []
    comment = ''

    # Read each line of the file
    for line in file_pointer.readlines():
        # Check if the line is not a comment
        if line[0] not in comment_indicators:
            line_data = []

            # Check if we're reading the full row or a subset of columns
            if row_len == -1:
                for d in line.split():
                    # Check if the data value is a NaN
                    if d == 'nan':
                        line_data.append(np.nan)
                    else:
                        line_data.append(float(d))
            else:
                for d in line.split()[:row_len]:
                    # Check if the data value is a NaN
                    if d == 'nan':
                        line_data.append(np.nan)
                    else:
                        line_data.append(float(d))

            data.append(line_data)
        else:
            # If the line is a comment, add it to the comment string
            comment += line
    
    # Close the file
    file_pointer.close()
    
    # Return the data array and the comment string
    return np.array(data), comment
