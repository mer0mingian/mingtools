# coding=utf-8
"""
code functions to export and read excitability data from/to csv-files.
"""
import numpy as np
a = numpy.asarray([ [1,2,3], [4,5,6], [7,8,9] ])


def write_excitability_to_file(filename, mdarray):
    """
    shoud touch filename and insert array mdarray rowwise.
    """
    np.savetxt(filename, mdarray, delimiter=",", format='%10.5f')


# with open(filename, 'w', 'a') as f:
#     f.write("\n".join(" ".join(map(str, x)) for x in (a,b)))



#def read_excitability_from_file(filename):
