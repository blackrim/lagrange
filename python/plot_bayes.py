import sys
import matplotlib.pyplot as plt


if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: python plot_bayes.py BAYES_OUT_ALL"
		sys.exit(0)
	infile = open(sys.argv[1],"r")
	ext = []
	labels = ["WNA->ENA","WNA->EU","WNA->EA","ENA->WNA","ENA->EU","ENA->EA","EU->WNA","EU->ENA","EU->EA","EA->WNA","EA->ENA","EA->EU"]
	params = {}
	for i in labels:
		params[i] = []
	for i in infile:
		spls = i.strip().split("\t")
		ext.append(float(spls[3]))
		start = 4
		for j in labels:
			params[j].append(float(spls[start]))
			start += 1
	plotparas = []
	for j in labels:
		plotparas.append(params[j])
	plt.boxplot(plotparas,positions=range(len(labels)))
	locs, labels = plt.xticks(range(len(labels)), labels)
	plt.show()
	infile.close()
