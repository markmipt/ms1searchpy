from sys import argv
import sys
sys.path.insert(0, '/home/mark/work/PycharmProjects/ms1searchpy/')
import main

inputfile = argv[1]#'/home/mark/work/PycharmProjects/ms1searchpy/tests/s6.mzML'
settings = main.settings(argv[2 ])

main.process_file(inputfile, settings)
print 'The search is finished.'