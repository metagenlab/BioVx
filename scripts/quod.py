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
import numpy as np


VERBOSITY = 0

class Wedge(object):
	def __init__(self, x, dx=0, y=0):
		self.x = x
		self.dx = dx
		self.y = y
	def draw(self, plt):
		(plt.xlim(), plt.ylim())
		if self.y > 0: ymin, ymax = 0.5, 1
		elif self.y < 0: ymin, ymax = 0, 0.5
		else: ymin, ymax = 0, 1
		x = max(plt.xlim()[0], min(plt.xlim()[1], self.x))
		plt.axvline(x=x, color='black', ymin=ymin, ymax=ymax)

		if self.dx:
			x = x + abs(self.dx)**.5 * (self.dx)/abs(self.dx)
			#if self.dx < 0: 
			if self.y == 0: y = 2
			else: y = self.y

			size = abs(self.dx)

			if self.dx < 0: marker = '<'
			else: marker = '>'

			plt.scatter((x,), (y,), marker=marker, color='black', s=25*size)

def hydro(gseq, window=19):
	#Copied with barely any modification from gblast3.py

	#sanitize the sequence, which may be a FASTA
	if gseq.startswith('>'):
		gseq = gseq[gseq.find('\n')+1:]
		gseq = re.sub('[^A-Z\-]', '', gseq)

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
			 'U':(0,0),
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
		return np.array(hydro)
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
	return np.array(newhydro)

try: 
	import huffman

	def entropy(gseq, window=19):
		if gseq.startswith('>'):
			gseq = gseq[gseq.find('\n')+1:]
			gseq = re.sub('[^A-Z\-]', '', gseq)
		seq=gseq.replace('-','')

		prev = 0
		midpt = (window+1)//2
		length = len(seq)
		entropy = []

		for i in range(length-window+1):
			h = huffman.Huffman(seq[i:i+window+1])
			entropy.append(h.compress())

		if len(seq) == len(gseq): return np.array(entropy)

		replace = re.finditer('(-+)',gseq)
		inserts = {}

		for i in replace: inserts.setdefault(i.start(),i.end()-i.start())

		first = False
		newentropy = []

		for x, h in enumerate(entropy):
			if x in inserts.keys() and first is False:
				first = True
				for y in range(inserts[x]):
					newentropy.append(matplotlib.numpy.nan)
				newcount = x + inserts[x]
				continue
			if first is False:
				newentropy.append(h)
				continue
			if first is True and newcount in inserts.keys():
				for y in range(inserts[newcount]):
					newentropy.append(matplotlib.numpy.nan)
				newcount += inserts[newcount]
				continue
			else:
				newentropy.append(h)
				newcount +=1
		return np.array(newentropy)
except ImportError: entropy = hydro

def hmmtop(sequence, silent=False):
	#Kevin's standard HMMTOP invocation

	sequence = sequence.replace('-', '')
	if not sequence.startswith('>'): sequence = '>seq\n' + sequence
	f = subprocess.Popen(['hmmtop', '-sf=FAS', '-is=pseudo', '-pi=spred', '-if=--'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	hmmtopout, err = f.communicate(input=sequence)
	print(err.strip(), file=sys.stderr)

	indices = map(int, re.findall('((?:[0-9]+\s+)+)$', hmmtopout)[0].split()[1:])
	pureseq = sequence[sequence.find('\n')+1:]
	pureseq = re.sub('[^A-Z\-]', '', pureseq)

	for i, c in enumerate(pureseq):
		if c in 'ACDEFGHIKLMNPQRSTVWY': pass
		else: 
			for j, n in enumerate(indices):
				if (n - 1) >= i: indices[j] += 1

	tmss = []
	for i in range(0, len(indices), 2): tmss.append(np.array(indices[i:i+2]))
	return tmss

def hydro_color(n):
	#make graphs less confusing (-K)
	if type(n) is int:
		r = n % 3
		if r == 0: return 'red'
		elif r == 1: return 'blue'
		elif r == 2: return 'green'
		#probably add back extra colors later?
	else: return n

def tms_color(n):
	#colors selected by Arturo
	if type(n) is int:
		r = n % 3
		if r == 0: return 'orange'
		elif r == 1: return 'cyan'
		elif r == 2: return 'darkgreen'
	else: return n

def what(sequences, labels=None, imgfmt='png', directory=None, filename=None, title=False, dpi=80, hide=True, viewer=None, bars=[], color='auto', offset=0, statistics=False, overwrite=False, manual_tms=None, wedges=[], ywedge=2, legend=False, window=19, silent=False, axisfont=None, tickfont=None, xticks=None, mode='hydropathy', width=None, height=None):
	#wedges: [(x1, dx1), (x2, dx2), ...]

	try: color = int(color)
	except ValueError: pass

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
	maxl = None
	hydropathies = []
	top = []

	for seq in sequences:
		if mode == 'hydropathy': hydropathies.append(hydro(seq, window))
		if mode == 'entropy': hydropathies.append(entropy(seq, window))


		if manual_tms: 
			tmss = []
			tms = manual_tms.pop(0)
			if tms[0]: tmss += (tms)
			else: tmss = [(tms, color) for tms in hmmtop(seq, silent=silent)]
			top.append(tmss)
		else:
			top.append([(tms, color) for tms in hmmtop(seq, silent=silent)])

		if minl is None: minl = len(hydropathies[-1])
		elif len(hydropathies[-1]) < minl: minl = len(hydropathies[-1])

		if maxl == None: maxl = len(hydropathies[-1])
		elif len(hydropathies[-1]) > maxl: maxl = len(hydropathies[-1])

	halen = len(hydropathies[0])

	#...except this, which may be an artifact of only worrying about 2 sequences...
	#X = range(len(hydropathies[0]))
	X = np.array(range(maxl))

	plt.figure()
	plt.axhline(y=0, color='black', linewidth=0.5)
	if mode == 'entropy': 
		plt.ylim(2, 5)
		plt.axhline(y=4.32193, color='black', linewidth=0.5)
	else: plt.ylim(-3, 3)
	#...and this one too...
	#plt.xlim(right=len(sequences[0]))
	plt.xlim(left=offset, right=maxl+offset+window)
	if xticks: plt.xticks(np.arange(offset, maxl+offset+window, xticks))
	if tickfont: plt.tick_params(labelsize=tickfont)
	if title: plt.suptitle(title)
	else:
		if len(labels) == 1: plt.suptitle(labels[0])
		elif len(labels) == 2: plt.suptitle('%s (red) and %s (blue)' % tuple(labels))
		elif len(labels) == 3: plt.suptitle('%s, %s, and %s' % tuple(labels))
		elif len(labels) > 3: plt.suptitle('%s, %s, %s, and %d more' % (tuple(labels[0:3]) + (len(labels) - 3,)))
	if axisfont:
		plt.xlabel('Residue number', fontsize=axisfont)
		if mode == 'entropy': plt.ylabel('Entropy', fontsize=axisfont)
		else: plt.ylabel('Hydropathy', fontsize=axisfont)
	else:
		plt.xlabel('Residue number')
		if mode == 'entropy': plt.ylabel('Entropy', fontsize=axisfont)
		else: plt.ylabel('Hydropathy', fontsize=axisfont)

	for i, seq in enumerate(sequences):
		hseq = hydropathies[i]

		if color == 'auto': c = i
		else: c = color
		plt.plot(X[:len(hseq)]+offset+window/2, hseq, linewidth=1, label=labels[i], color=hydro_color(c))
		for tms in top[i]: 
			if tms[1] == 'auto': c = i
			else: c = tms[1]
			plt.axvspan(tms[0][0]+offset, tms[0][1]+offset, facecolor=tms_color(c), alpha=0.25)

		#...and last but greatest, this, which eliminates crucial logic from gblast3.py
		#I can only assume it makes sure TMSs land on the right spots in alignments with lots of -'s and Xs
		#Mostly, I don't quite know how self.queryhmg or self.tcdbhmg return their results
		#for tms in top[i]: plt.axvspan(tms[0], tms[1], facecolor=color(i, color), alpha=0.6/len(hydropathies))

	#if bars:
	#	for x in bars: plt.axvline(x=x, color='red', ymin=0.5, ymax=1)
	#if wedges:
	#	for wedge in wedges:
	#		wedge[0] += abs(wedge[1])**.5 * (wedge[1])/abs(wedge[1])
	#		if wedge[1] < 0: 
	#			wedge[1] *= -1
	#			marker = '<'
	#		else:
	#			marker = '>'
	#		plt.scatter((wedge[0],), (ywedge,), marker=marker, color='black', s=25*wedge[1])
	for w in wedges: w.draw(plt)
	
	plt.legend().set_visible(legend)
	#fig = plt.gcf()
	#fig.set_size_inches(15, 3)

	if width == None: width = (0.0265)*maxl if maxl > 600 else 15
	if height == None: height = 5.5
	plt.gcf().set_size_inches(width, height)

	plt.savefig(filename, dpi=dpi, format=imgfmt, bbox_inches='tight', pad_inches=0.003)
	if VERBOSITY != VERBOSITY: 
		if len(labels) == 1: print('%s: %s' % (filename, labels[0]))
		elif len(labels) == 2: print('%s: %s, %s' % (filename, labels[0], labels[1]))
		elif len(labels) == 3: print('%s: %s, %s, %s' % (filename, labels[0], labels[1], labels[2]))
		elif len(labels) > 3: print('%s: %s, %s, %s, and %d others' % (filename, labels[0], labels[1], labels[2], len(labels)-3))
		else: print(filename)
	if not hide and (imgfmt != 'eps' and imgfmt != 'tif'):
		if viewer: IMAGE_VIEWER = viewer
		else:
			if sys.platform.startswith('linux'): IMAGE_VIEWER = 'xdg-open'
			elif sys.platform == 'darwin': IMAGE_VIEWER = 'open'
			else: print('[WARNING]: Unknown platform; Email khendarg@ucsd.edu with your OS and default generic file-opening program', file=sys.stderr)
		os.system('%s %s' % (IMAGE_VIEWER, filename))
	plt.clf()
	plt.close()
	#fig.clear()

def parse_csranges(colranges):
	#assumes colranges is already in ['x-y:color', 'i-j:color,u-v:color'] format
	if not colranges: return []
	out = []
	colors = []
	for rangeset in colranges:
		out.append([])
		colors.append([])
		if rangeset.strip().lower().startswith('skip'):
			sspan = re.split('\s*:\s*', rangeset)
			if len(sspan) < 2: color = 'auto'
			else:
				try: color = int(sspan[1])
				except ValueError: color = sspan[1]
			out[-1] = None, color
			continue
		for span in re.split('\s*,\s*', rangeset):
			sspan = re.split('\s*:\s*', span)
			if len(sspan) < 2: color = 'auto'
			else:
				try: color = int(sspan[1])
				except ValueError: color = sspan[1]
			out[-1].append(([int(x) for x in re.split('\s*-\s*', sspan[0])], color))

	return out
			

def parse_wranges(wranges):
	#assumes wranges is already in ['x1,dx1,y1', x2,dx2,y2'] format
	if not wranges: return []
	out = []
	for wedge in wranges:
		props = [float(x) for x in re.split('\s*,\s*', wedge)]

		x, dx, y = 0, 0, 0
		if len(props) >= 1: x = props[0]
		if len(props) >= 2: dx = props[1]
		if len(props) >= 3: y = props[2]
		out.append(Wedge(x, dx, y))

	return out

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='QUOD, the Questionable Utility of Doom: Turns sequences into WHAT plots, saves them in a convenient directory (like a criminally empty %s), and opens them' % os.getcwd())

	parser.add_argument('inp', metavar='input', nargs='+', help='Input files. Use - for stdin, also known as pasting the sequences in separated by ^Ds. Without -s or -f or with both -s and -d, inputs are automatically interpreted as either sequences or filenames.')
	parser.add_argument('-a', metavar='viewer', default=None, help='Viewer to be used for opening graphs')
	parser.add_argument('-b', '--bars', nargs='+', type=int, help='Draws vertical bars at these positions')
	parser.add_argument('-B', '--box', nargs='+', type=int, help='Draws boxes at these positions')
	parser.add_argument('-c', '--color', metavar='color', default='auto', help='Colors all curves')
	parser.add_argument('-d', metavar='directory', default=None, help='Directory to store graphs in (recommended only with autogenerated filenames)')
	parser.add_argument('-e', action='store_true', help='switch to entropy graphing mode')
	parser.add_argument('-f', action='store_true', help='Force all inputs to be interpreted as filenames (may be changed)')
	#TODO
	parser.add_argument('-i', default=None, help='File containing drawing instructions')
	parser.add_argument('-l', metavar='graph_title', help='Label graph with a specific title')
	parser.add_argument('-m', '--manual-tms', metavar='TMSs', nargs='+', help='Use these comma-separated ranges as TMSs instead of HMMTOP output. Use spaces to separate ranges for other sequences and \'skip\' to skip sequences (i.e. letting HMMTOP assign TMSs. Use colons to specify colors for individual TMSs. Example: -m 1-20:orange paints residues between 1 and 20 inclusive orange. -m 30-50,55-75:cyan paints residues between 30 and 50 inclusive the default color and residues between 55 and 75 inclusive cyan.')
	parser.add_argument('-o', metavar='filename', default=None, help='Filename of graph, relative to the argument of -d if present and as normally interpreted otherwise')
	parser.add_argument('-q', action='store_true', help='"quiet" mode, disables automatic opening of graphs')
	parser.add_argument('-r', metavar='dpi', default=80, type=int, help='Resolution of graph in dpi. The default is 80dpi, which is suitable for viewing on monitors. Draft-quality images should be about 300dpi, and publication-quality images need to be 600 or 1200dpi depending on the journal.')
	parser.add_argument('-s', action='store_true', help='Force inputs to be interpreted as sequences (may be changed)')
	parser.add_argument('-t', metavar='format', default='png', help='Format of graph; the default is png. Also accepted are eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, and tiff.')
	parser.add_argument('-v', action='store_true', help='Verbose output. Enables warnings and generally makes things messier')
	parser.add_argument('-W', '--walls', metavar='x(,dx(,y))', nargs='+', help='Shorthand to specifying both --wedges and --bars for the same x-values. Cumulative with both')
	parser.add_argument('-w', '--wedges', metavar='wedgex,wedgedx', nargs='+', help='Draw dx-long wedges starting at x. Negative dx values result in left-pointing wedges, and positive dx values result in right-pointing wedges.')
	parser.add_argument('--axis-font', metavar='size', type=int, help='Axis label size')
	parser.add_argument('--height', metavar='height', type=float, help='Plot height in inches {default:5.5}')
	parser.add_argument('--offset', metavar='init_resi', default=0, type=int, help='Sets starting x-value')
	parser.add_argument('--statistics', action='store_true', help='Display some useful statistics on the results')
	parser.add_argument('--tick-font', metavar='size', type=int, help='Tick label size')
	parser.add_argument('--width', metavar='width', type=float, help='Plot width in inches {default:dynamic}')
	parser.add_argument('--xticks', metavar='interval', type=float, help='Spacing between x ticks')
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

	if args.e: mode = 'entropy'
	else: mode = 'hydropathy'

	#wedges, bars = [], []
	#
	#if args.walls: 
	#	wedges += parse_wranges(args.walls)
	#	bars += [x[0] for x in wedges]
	#if args.wedges: wedges += parse_wranges(args.wedges)
	#if args.bars: bars += [int(x) for x in args.bars]

	wedges = []
	if args.walls:
		wedges += parse_wranges(args.walls)
	if args.wedges:
		wedges += parse_wranges(args.wedges)
	if args.bars:
		wedges += [Wedge(int(x), 0, 0) for x in args.bars]

	what(sequences, labels=labels, imgfmt=args.t, directory=args.d, filename=args.o, title=args.l, dpi=args.r, hide=args.q, viewer=args.a, overwrite=True, color=args.color, statistics=args.statistics, offset=args.offset, manual_tms=parse_csranges(args.manual_tms), wedges=wedges, xticks=args.xticks, axisfont=args.axis_font, tickfont=args.tick_font, mode=mode, width=args.width, height=args.height)
