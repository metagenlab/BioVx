#!/usr/bin/env python2
#QUOD: Questionable Utility Of Doom
#most code copied from gblast3.py (which was written by Vamsee Reddy and Vasu Pranav Sai Iddamsetty) except where noted
#-Kevin Hendargo
from __future__ import print_function, division

from Bio import SeqIO
import argparse
import tempfile, os, sys, re
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

VERBOSITY = 0

def hydro(gseq, window=19):
	#Copied with barely any modification from gblast3.py
	seq=gseq.replace('-','')
	prev = 0
	index = {'G':(-0.400,0.48),
			 'I':(4.500,1.38),
			 'S':(-0.800,-0.18),
			 'Q':(-3.500,-0.85),
			 'E':(-3.500,-0.74),
			 'A':(1.800,0.62),
			 'M':(1.900,0.64),
			 'T':(-0.700,-0.05),
			 'Y':(-1.300,0.26),
			 'H':(-3.200,-0.4),
			 'V':(4.200,1.08),
			 'F':(2.800,1.19),
			 'C':(2.500,0.29),
			 'W':(-0.900,0.81),
			 'K':(-3.900,-1.5),
			 'L':(3.800,1.06),
			 'P':(-1.600,0.12),
			 'N':(-3.500,-0.78),
			 'D':(-3.500,-0.90),
			 'R':(-4.500,-2.53),
			 'X':(0,0),
			 'B':(-3.500,-0.84),
			 'J':(-3.500,-0.80),
			 'Z':(4.150,1.22)
			}
			#B, J, Z: arithmetic means of [D,N], [I,L], and [E,Q] respectively (-Kevin)
	midpt = (window+1)//2
	length = len(seq)
	hydro = []
	for i in range(length-window+1):
		total = 0
		for j in range(window):
			total += index[seq[i+j]][0]
		total = total/window
		hydro.append(total)
	if len(seq) == len(gseq):
		return hydro
	replace = re.finditer('(-+)',gseq)
	inserts = {}
	for i in replace:
		inserts.setdefault(i.start(),i.end()-i.start())
	first = False
	newhydro = []
	for x, h in enumerate(hydro):
		if x in inserts.keys() and first is False:
			first = True
			for y in range(inserts[x]):
				newhydro.append(matplotlib.numpy.nan)
			newcount = x + inserts[x]
			continue
		if first is False:
			newhydro.append(h)
			continue
		if first is True and newcount in inserts.keys():
			for y in range(inserts[newcount]):
				newhydro.append(matplotlib.numpy.nan)
			newcount += inserts[newcount]
			continue
		else:
			newhydro.append(h)
			newcount +=1
	return newhydro

def hmmtop(sequence):
	#Kevin's standard HMMTOP invocation
	f = tempfile.NamedTemporaryFile(delete=False)
	try:
		f.write('>tmpsequence\n' + sequence)
		f.close()
		hmmtopout = subprocess.check_output(['hmmtop', '-sf=FAS', '-is=pseudo', '-pi=spred', '-if=' + f.name], stderr=open(os.devnull, 'w'))
	finally: os.remove(f.name)

	out = hmmtopout.split()
	ntmss = int(out[4])
	if not ntmss: return []

	tmss = []
	for i in range(ntmss): tmss.append(map(int, out[5 + 2*i:7+2*i]))
	return tmss

def hydro_color(n):
	#make graphs less confusing (-K)
	r = n % 3
	if r == 0: return 'red'
	elif r == 1: return 'blue'
	elif r == 2: return 'green'
	#probably add back extra colors later?

def tms_color(n):
	#colors selected by Arturo
	r = n % 2
	if r == 0: return 'orange'
	elif r == 1: return 'cyan'

def what(sequences, labels=None, imgfmt='png', directory=None, filename=None, title=False, dpi=80, hide=True, viewer=None, overwrite=False, window=19):
	#generalized from gblast3.py but mostly the same...
	if not labels:
		labels = []
		for i in range(len(sequences)): labels.append('Sequence %d' % i)

	#make directories look nice
	if directory and not directory.endswith('/'): directory += '/'

	#autogenerate a filename if none was specified
	if filename == None:
		if directory: filename = directory
		else: filename = ''
		filename += re.sub('[^A-Za-z_\.\-0-9]', '_', labels[0].replace(' ', '_'))
		try:
			for l in labels[1:4]: filename += '_' + re.sub('[^A-Za-z_\.\-0-9]', '_', re.sub('[\s/]', '_', l))
		except IndexError: pass
		if len(labels) > 3: filename += '_etc'

		filename += '.' + imgfmt
	#set filename if directory and filename were specified
	elif filename and directory:
		filename = directory + filename
	#else, just use the filename as is

	if not overwrite and os.path.isfile(filename): return

	minl = None
	#maxl = None
	maxl = 0
	hydropathies = []
	top = []

	for seq in sequences:

		top.append(hmmtop(seq))
		
		hydropathies.append(hydro(seq, window=window))

		if minl == None: minl = len(hydropathies[-1])
		elif len(hydropathies[-1]) < minl: minl = len(hydropathies[-1])

		#if maxl == None: maxl = len(hydropathies[-1])
		#elif len(hydropathies[-1]) > maxl: maxl = len(hydropathies[-1])
		if len(seq) > maxl: maxl = len(seq)

	halen = len(hydropathies[0])

	#...except this, which may be an artifact of only worrying about 2 sequences...
	#X = range(len(hydropathies[0]))
	X = range(maxl)

	plt.figure()
	plt.axhline(y=0, color='black')
	plt.ylim(-3, 3)
	#...and this one too...
	#plt.xlim(right=len(sequences[0]))
	plt.xlim(right=maxl)

	if title: plt.suptitle(title)

	else:
		if len(labels) == 1: plt.suptitle(labels[0])
		elif len(labels) == 2: plt.suptitle('%s (red) and %s (blue)' % tuple(labels))
		elif len(labels) == 3: plt.suptitle('%s, %s, and %s' % tuple(labels))
		elif len(labels) > 3: plt.suptitle('%s, %s, %s, and %d more' % (tuple(labels[0:3]) + (len(labels) - 3,)))

	plt.xlabel('Residue #')
	plt.ylabel('Hydro')
	plt.legend(loc='lower right')

	for i, seq in enumerate(sequences):
		hseq = hydropathies[i]

		plt.plot(X[:len(hseq)], hseq, linewidth=1, label=labels[i], color=hydro_color(i))

		#...and last but greatest, this, which eliminates crucial logic from gblast3.py
		#I can only assume it makes sure TMSs land on the right spots in alignments with lots of -'s and Xs
		#Mostly, I don't quite know how self.queryhmg or self.tcdbhmg return their results
		#for tms in top[i]: plt.axvspan(tms[0], tms[1], facecolor=color(i), alpha=0.6/len(hydropathies))
		for tms in top[i]: 
			plt.axvspan(tms[0]-window/2, tms[1]-window/2, facecolor=tms_color(i), alpha=0.25)

	fig = plt.gcf()
	fig.set_size_inches(15, 3)
	#f = tempfile.NamedTemporaryFile()
	#f.close()
	plt.savefig(filename, dpi=dpi, format=imgfmt, bbox_inches='tight', pad_inches=0.003)
	#if VERBOSITY != VERBOSITY: 
	#	if len(labels) == 1: print('%s: %s' % (filename, labels[1]))
	#	elif len(labels) == 2: print('%s: %s, %s' % (filename, labels[0], labels[1]))
	#	elif len(labels) == 3: print('%s: %s, %s, %s' % (filename, labels[0], labels[1], labels[2]))
	#	elif len(labels) > 3: print('%s: %s, %s, %s, and %d others' % (filename, labels[0], labels[1], labels[2], len(labels)-3))
	#        else: print(filename)
	if not hide and (imgfmt != 'eps' and imgfmt != 'tif'):
		if viewer: IMAGE_VIEWER = viewer
		else:
			if sys.platform.startswith('linux'): IMAGE_VIEWER = 'xdg-open'
			elif sys.platform == 'darwin': IMAGE_VIEWER = 'open'
			else: print('[WARNING]: Unknown platform; Email khendarg@ucsd.edu with your OS and default generic file-opening program', file=sys.stderr)
		os.system('%s %s' % (IMAGE_VIEWER, filename))
	plt.clf()
	plt.close()
	fig.clear()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='QUOD, the Questionable Utility of Doom: Turns sequences into WHAT plots, saves them in a convenient directory (like a criminally empty %s), and opens them' % os.getcwd())

	parser.add_argument('inp', metavar='input', nargs='+', help='Input files. Use - for stdin, also known as pasting the sequences in separated by ^Ds. Without -s or -f or with both -s and -d, inputs are automatically interpreted as either sequences or filenames.')
	parser.add_argument('-f', action='store_true', help='Force all inputs to be interpreted as filenames (may be changed)')
	parser.add_argument('-s', action='store_true', help='Force inputs to be interpreted as sequences (may be changed)')
	parser.add_argument('-v', action='store_true', help='Verbose output. Enables warnings and generally makes things messier')
	parser.add_argument('-d', metavar='directory', default=None, help='Directory to store graphs in (recommended only with autogenerated filenames)')
	parser.add_argument('-o', metavar='filename', default=None, help='Filename of graph, relative to the argument of -d if present and as normally interpreted otherwise')
	parser.add_argument('-q', action='store_true', help='"quiet" mode, disables automatic opening of graphs')
	parser.add_argument('-a', metavar='viewer', default=None, help='Viewer to be used for opening graphs')
	parser.add_argument('-t', metavar='format', default='png', help='Format of graph; the default is png. Also accepted are eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, and tiff.')
	parser.add_argument('-r', metavar='dpi', default=80, type=int, help='Resolution of graph in dpi. The default is 80dpi, which is suitable for viewing on monitors. Draft-quality images should be about 300dpi, and publication-quality images need to be 600 or 1200dpi depending on the journal.')
	parser.add_argument('-l', metavar='graph_title', help='Label graph with a specific title')
	args = parser.parse_args()

	#warnings look messy
	if args.v: VERBOSITY = 1
	if not VERBOSITY:
		import warnings
		warnings.filterwarnings('ignore', '.')

	sequences = []
	labels = []
	n = 0
	for inp in args.inp:
		#read stuff from stdin
		if inp == '-':
			seq = ''
			try:
				while 1: seq += raw_input()
			except EOFError: 
				sequences.append(seq)
				labels.append('Sequence %d' % n)
		elif not(args.f ^ args.s):
			#does a dumb autodetect if a user uses both or neither of -f and -s 
			if re.findall('[^A-Z\-]', inp): 
				for seq in SeqIO.parse(inp, 'fasta'): 
					sequences.append(str(seq.seq))
					labels.append(seq.name)
			else: 
				sequences.append(inp)
				labels.append('Sequence %d' % n)

		elif args.f:
			for seq in SeqIO.parse(inp, 'fasta'): 
				sequences.append(str(seq.seq))
				labels.append(seq.name)

		elif args.s: 
			sequences.append(inp)
			labels.append('Sequence %d' % n)
		n += 1
	what(sequences, labels=labels, imgfmt=args.t, directory=args.d, filename=args.o, title=args.l, dpi=args.r, hide=args.q, viewer=args.a, overwrite=True)
