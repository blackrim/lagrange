import sys

def range_string(num_ranges,range_tuple):
	range_list = ["0"]*int(num_ranges)
	for i in range_tuple:
		range_list[int(i)] = "1"
	range_str = "".join(range_list)
	return range_str

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "python convert_configurator_lagrange_cpp.py infile.py outbasename"
		sys.exit(0)
	
	infile = open(sys.argv[1],"r")
	read = False
	data = ""
	for i in infile:
		if "begin data" in i: 
			read = True
		if read == True:
			data += i
			if "end data" in i:
				read = False
				break
	infile.close()
	d = eval(data)
	print "data successfully read from "+sys.argv[1]
	#print d
	outbasename = sys.argv[2]
	print "creating lagrange cpp files"
	#print the trees to a file
	treefile = open(outbasename+".tre","w")
	for i in d['newick_trees']:
		treefile.write(i['newick']+"\n")
	treefile.close()
	#print the data to a file
	datafile = open(outbasename+".data","w")
	datafile.write(str(len(d['taxon_range_data']))+"\t"+str(len(d['area_labels']))+"\n")
	for i in d['taxon_range_data']:
		datafile.write(i+"\t"+range_string(len(d['area_labels']),d['taxon_range_data'][i])+"\n")
	datafile.close()
	#print the main run file
	mainfile = open(outbasename+".lg","w")
	mainfile.write("treefile = "+outbasename+".tre\n")
	mainfile.write("datafile = "+outbasename+".data\n")
	mainfile.write("ratematrix = "+outbasename+".rm\n")
	mainfile.write("areanames = "+str(" ".join(d['area_labels']))+"\n")
	mainfile.write("splits\n")
	mainfile.write("ancstate = _all_\n")
	mainfile.write("periods = ")
	curtime = 0
	for i in d['dispersal_durations']:
		mainfile.write(str(float(i)+curtime))
		curtime += float(i)
		if i != d['dispersal_durations'][-1]:
			mainfile.write(" ")
	mainfile.write("\n")
	mainfile.write("includedists = "+str(d['max_range_size'])+"\n")
	mainfile.write("includedists = ")
	for i in d['ranges']:
		mainfile.write(range_string(len(d['area_labels']),i))
		if i != d['ranges'][-1]:
			mainfile.write(" ")
	mainfile.write("\n")
	mainfile.write("excludedists = ")
	for i in d['excluded_ranges']:
		mainfile.write(range_string(len(d['area_labels']),i))
		if i != d['excluded_ranges'][-1]:
			mainfile.write(" ")
	mainfile.write("\n")
	mainfile.write("estimate_dispersal_mask\n")
	mainfile.close()
	#print the rate matrix file
	ratefile = open(outbasename+".rm","w")
	for i in d['area_dispersal']:
		for j in i:
			sj = j
			for k in range(len(sj)):
				sj[k] = str(sj[k])
			ratefile.write(" ".join(sj)+"\n")
		ratefile.write("\n")
	ratefile.close()
	print "creation finished"

	print "-----"

	print "you should be able to run the analyses now with something like"
	print "lagrange "+str(outbasename)+".lg"
