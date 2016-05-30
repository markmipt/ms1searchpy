from sys import argv
import main

inputfile = argv[1]
settings = main.settings(argv[2 ])

main.process_file(inputfile, settings)
print 'The search is finished.'