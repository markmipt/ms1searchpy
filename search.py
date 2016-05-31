from sys import argv
import main, utils

inputfile = argv[1]
settings = utils.settings(argv[2 ])

main.process_file(inputfile, settings)
print 'The search is finished.'