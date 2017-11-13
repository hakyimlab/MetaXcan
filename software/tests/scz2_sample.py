import pandas
import numpy

expected_snp = pandas.Series(["rs940550", "rs6650104", "rs6594028", "rs9701055", "rs7417504", "rs12082473", "rs3094315", "rs3131971", "rs61770173", ], dtype=numpy.str)
expected_effect = pandas.Series(["C", "T", "T", "A", "T", "A", "A", "T", "A"], dtype=numpy.str)
expected_non_effect = pandas.Series(["G", "C", "C", "T", "C", "G", "G", "C", "C"], dtype=numpy.str)
expected_chromosome = pandas.Series(["chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr22", "chr22"], dtype=numpy.str)
expected_position = pandas.Series([729679, 731718, 734349, 736289, 751756, 752566, 752721, 752894, 753405])
expected_se = pandas.Series( [0.0173, 0.0198, 0.02, 0.0193, 0.0164, 0.0149, 0.0146, 0.015, 0.0159], dtype=numpy.float32)
expected_p = pandas.Series( [0.2083, 0.3298, 0.3055, 0.5132, 0.8431, 0.7870, 0.8229, 0.5065, 0.8181], dtype=numpy.float32)

expected_beta = pandas.Series( [-0.0217038334437866, 0.0193025022544974, 0.0204984635773248, -0.0125990355765152, -0.00319509889654086, -0.00399798128729837, 0.00330453400830047, 0.0099998345783334, -0.00369682484428976], dtype=numpy.float32)
expected_beta_fix = pandas.Series( [-0.02170383,  0.0193025 ,  0.02049846, -0.01259904, -0.0031951, -0.00399798,  0.00330453,  0.00999983,  2.22032681], dtype=numpy.float32)

expected_zscore_1 = pandas.Series( [-1.254557, 0.974874, 1.024923, -0.652800, -0.194823, -0.268321, 0.226338, 0.666656, -0.232505], dtype=numpy.float32)
expected_zscore_2 = pandas.Series( [-1.258254,  0.974517,  1.02471 , -0.653863, -0.19793 , -0.270208, 0.223817,  0.664297, -0.229989], dtype=numpy.float32)
expected_zscore_fix = pandas.Series( [ -1.25825381,   0.97451681,   1.02471006,  -0.65386313, -0.19792978,  -0.2702083 ,   0.22381662,   0.6642974 ,  11.52388382], dtype=numpy.float32)
