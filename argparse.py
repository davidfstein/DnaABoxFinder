import sys, getopt

def parseargs(argv):
	genome_file_path = ''
	window_length = 0
	num_kmers = 0
	dnaa_box_length = 0
	hamming_distance = 0
	
	try:
		opts, args = getopt.getopt(argv, "p:w:n:d:h:", ["path=", 
													    "window-length=", 
													    "num-kmers=", 
													    "dnaa-box-length", 
													    "hamming-distance="])
	except getopt.GetoptError:
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-p' or opt == '--path':
			genome_file_path = arg
		if opt == '-w' or opt == '--window-length':
			window_length = int(arg)
		if opt == '-n' or opt == '--num-kmers':
			num_kmers = int(arg)
		if opt == '-d' or opt == '--dnaa-box-length':
			dnaa_box_length = int(arg)
		if opt == '-h' or opt == '--hamming-distance':
			hamming_distance = int(arg)
	return genome_file_path, window_length, num_kmers, dnaa_box_length, hamming_distance
