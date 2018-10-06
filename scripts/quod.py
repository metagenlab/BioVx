#!/usr/bin/env python2
#QUOD2: Questionable Utility of Doom 2
#most code copied from gblast3.py (which was written by Vamsee Reddy and Vasu Pranav Sai Iddamsetty) except where noted
#-Kevin Hendargo
from __future__ import division, print_function, division
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.patches as patches
import os, subprocess, re, sys
import argparse
import numpy as np

from Bio import SeqIO

def info(*things): print('[INFO]:', *things, file=sys.stderr)

VERBOSITY = 0
MISTAKEFUDGE = 600
HYDROPATHY = {'G':-0.400, \
	 'I':4.500, \
	 'S':-0.800, \
	 'Q':-3.500, \
	 'E':-3.500, \
	 'A':1.800, \
	 'M':1.900, \
	 'T':-0.700, \
	 'Y':-1.300, \
	 'H':-3.200, \
	 'V':4.200, \
	 'F':2.800, \
	 'C':2.500, \
	 'W':-0.900, \
	 'K':-3.900, \
	 'L':3.800, \
	 'P':-1.600, \
	 'N':-3.500, \
	 'D':-3.500, \
	 'R':-4.500, \
	 'U':0, \
	 'B':-3.500, \
	 'J':-3.500, \
	 'Z':4.150 \
}
AMPHIPATHICITY = {'G':0.48, \
	 'I':1.38, \
	 'S':-0.18, \
	 'Q':-0.85, \
	 'E':-0.74, \
	 'A':0.62, \
	 'M':0.64, \
	 'T':-0.05, \
	 'Y':0.26, \
	 'H':-0.4, \
	 'V':1.08, \
	 'F':1.19, \
	 'C':0.29, \
	 'W':0.81, \
	 'K':-1.5, \
	 'L':1.06, \
	 'P':0.12, \
	 'N':-0.78, \
	 'D':-0.90, \
	 'R':-2.53, \
	 'U':0, \
	 'B':-0.84, \
	 'J':-0.80, \
	 'Z':1.22 \
}
CHARGE = {'G':0.0, \
	 'I':0.0, \
	 'S':0.0, \
	 'Q':0.0, \
	 'E':-1.0, \
	 'A':0.0, \
	 'M':0.0, \
	 'T':0.0, \
	 'Y':+0.001, \
	 'H':+0.1, \
	 'V':0.0, \
	 'F':0.0, \
	 'C':+0.1, \
	 'W':0.0, \
	 'K':+1.0, \
	 'L':0.0, \
	 'P':0.0, \
	 'N':0.0, \
	 'D':-1.0, \
	 'R':+1.0, \
	 'U':0, \
	 'B':-0.5, \
	 'J':-0.5, \
	 'Z':0.0 \
}

def isspans(token):
	''' checks if a token is a span/range '''
	if re.match(r'^[0-9]+(?:-[0-9]+)(?:,[0-9]+(?:-[0-9]+))*$', token): return True
	else: return False

def parse_spans(token):
	''' parses a token as a set of spans '''
	spans = []
	for span in token.split(','):
		indices = span.split('-')
		if len(indices) == 0: continue
		elif len(indices) == 1: spans.append([int(indices)] * 2)
		elif len(indices) >= 2: spans.append([int(index) for index in indices[:2]])
	return spans

def isid(token):
	''' checks if a token is an ID '''
	if re.match('^\+[0-9]+$', token): return True
	else: return False

def parse_id(token):
	''' parses an ID '''
	if isid(token): return int(token)
	else: return None

def isint(token):
	''' checks if something is an int '''
	try:
		int(token)
		return True
	except ValueError: return False

def isfloat(token):
	''' checks if something is a float '''
	try:
		float(token)
		return True
	except ValueError: return False

def warn(*args): 
	''' prints a warning '''
	print('[WARNING]:', *args, file=sys.stderr)

def sanitize(s): 
	''' removes forslashes to generate UNIX-compatible filenames '''
	return re.sub('/', '', s)

def roundto(x, n=5):
	''' round to the next lowest multiple of n '''
	return (x//n)*n

def overlap(span1, span2):
	''' checks if two spans overlap '''
	if span1[0] <= span2[0] <= span1[-1]: return True
	elif span1[0] <= span2[1] <= span1[-1]: return True
	elif span2[0] <= span1[0] <= span2[-1]: return True
	elif span2[0] <= span1[1] <= span2[-1]: return True
	else: return False

def union(span1, span2):
	''' produces unions of contiguous/intersecting spans '''
	if not overlap(span1, span2): raise NotImplementedError
	else: return [min(span1[0], span2[0]), max(span1[-1], span2[-1])]

def hex2tuple(s):
	''' turns hexes into tuples '''
	if s.startswith('#'): s = s[1:]
	elif s.startswith('0x'): s = s[2:]

	if len(s) == 3: s = [2*s[i] for i in range(len(s))]
	if len(s) == 1: s *= 6
	
	l = [int(s[i:i+2], 16)/255. for i in range(0, len(s), 2)]
	return tuple(l)

class Plot(object):
	''' a plot, as in a pair of axes and an area '''
	def __init__(self, fig=None, canvas=None): 
		''' constructor
		fig: Figure object
		canvas: FigureCanvas object
		'''
		self.fig = Figure() if fig is None else fig
		self.canvas = FigureCanvas(self.fig) if fig is None else canvas

		self.ax = self.fig.add_subplot(111)
		self.width, self.height = None, None
		self.x_offset = 0
		self.xlim = [0, 20]
		self.ylim = [-3, 3]
		self.xticks, self.yticks = None, 1

		self.axeslabels = ['X', 'Y']
		self.titles = []

		self.grid = False

		self.entities = []

	def add(self, entity): self.entities.append(entity)
	def pop(self, i): self.entities.pop(i)
	def remove(self, entity): self.entities.remove(entity)
	def __iter__(self): return iter(self.entities)
	def __getitem__(self, i): return self.entities[i]
	def __getslice__(self, s): return self.entities[s]

	def render(self):
		''' configure plot appropriately '''

		for e in self.entities: 
			try: e.draw(self)
			except AttributeError:
				print(e)
				raise AttributeError

		self.ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
		self.ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])
		xlim = self.ax.get_xlim()
		ylim = self.ax.get_ylim()

		if self.x_offset:
			self.ax.set_xlim([x + self.x_offset for x in xlim])

		xlim = self.ax.get_xlim()
		ylim = self.ax.get_ylim()

		maxl = xlim[1] - xlim[0]
		if self.xticks is None: self.xticks = roundto(maxl // 5, 5)

		self.ax.set_xticks(np.arange(xlim[0], xlim[1]+1, max(1, self.xticks)))
		self.ax.set_yticks(np.arange(ylim[0], ylim[1]+self.yticks, max(0.01, self.yticks)))

		#TODO: allow centimeters or something
		#if self.width is None: self.width = (0.0265) * (maxl) if maxl > 600 else 15.
		if self.width is None: self.width = (0.0265) * (maxl) if maxl > 200 else 15.
		if self.height is None: self.height = 5.5
		self.fig.set_figwidth(self.width)
		self.fig.set_figheight(self.height)
		self.fig.set_tight_layout(True)
		self.ax.set_xlabel(self.axeslabels[0])
		self.ax.set_ylabel(self.axeslabels[1])

		#only when plotting hydropathy! FIXME
		self.ax.axhline(y=0, color='black', linewidth=0.5)

		title = ''
		for t in self.titles: 
			if t is not None: title += '{}, '.format(t)
		if len(title) > 40: title = title[:40] + '.....'

		self.ax.set_title(title[:-2])

		if self.grid: self.ax.grid('on')

	def save(self, filename, dpi=80, format='png'):
		''' save plot to disk '''
		self.fig.savefig(filename, dpi=dpi, format=format)

class Entity(object):
	''' something that can be drawn on a plot '''
	def __init__(self): pass

	def draw(self, plot): 
		''' called when rendering plots '''
		raise NotImplementedError

class Curve(Entity):
	''' a curve '''
	def __init__(self, X=[], Y=[], style='auto'):
		''' constructor
		X: iterable of floats (default:[])
		Y: iterable of floats (default:[])
		style: int/str, QUOD/Matplotlib-compatible color specifier (default:'auto')
		'''
		Entity.__init__(self)
		self.X = X
		self.Y = Y
		self.style = style

	def set_color(self, style=None):
		''' set Curve color. If style is an int, automatically select from three '''
		try: style = int(style)
		except ValueError: pass
		except TypeError: pass

		if style is None: style = self.style

		if type(style) is int: 
			if (style % 3 == 0): self.style = 'red'
			if (style % 3 == 1): self.style = 'blue'
			if (style % 3 == 2): self.style = 'green'
		else: self.style = style

	def draw(self, plot): 
		''' called when rendering plots '''
		self.set_color()
		#plot.xlim[0] = min(plot.xlim[0], self.X[0])
		#plot.xlim[1] = max(plot.xlim[1], self.X[-1])
		if len(self.Y) and not len(self.X): 
			plot.ax.plot(self.Y, self.style)
		else: plot.ax.plot(self.X, self.Y, self.style)

	def convolve(self, other, **kwargs):
		''' convolve with another curve '''
		out = Curve(**kwargs)
		out.X = np.copy(self.X[:(abs(len(self.X) - len(other.X)) + 1)])
		out.Y = np.zeros(len(out.X))

		for i, y_in in enumerate(self.Y[:len(out.Y)]):
			for j, y_filt in enumerate(other.Y): 
				try: out.Y[i+j] += y_in * y_filt
				except IndexError: pass
		out.X -= len(other.Y)//2
		return out

class Vspans(Entity):
	''' things with a width and infinite height '''
	def __init__(self, spans=[], style='orange', alpha=None):
		''' constructor
		spans: iterable of iterables. (default:[])
		style: color specifier (default:'orange')
		alpha: transparency/opacity (default:None)
		'''
		Entity.__init__(self)
		self.spans = spans
		self.style = style
		self.alpha = alpha

	def get_alpha(self):
		''' get transparency '''
		if type(self.style) is tuple and len(self.style) == 4: return self.style[3]
		elif self.alpha is not None: return self.alpha
		else: return 0.25

	def set_color(self, style=None):
		''' set color or choose a default '''
		try: style = int(style)
		except ValueError: pass
		except TypeError: pass

		if style is None: style = self.style

		if type(style) is int: 
			if (style % 3 == 0): self.style = 'orange'
			if (style % 3 == 1): self.style = 'cyan'
			if (style % 3 == 2): self.style = 'purple'
		else: self.style = style

	def draw(self, plot):
		''' called when rendering plots '''
		self.set_color()
		for span in self.spans:
			plot.ax.axvspan(span[0], span[1], facecolor=self.style, alpha=self.get_alpha())

class Region(Vspans):
	''' little boxes with widths, y-positions, and text but no height '''
	def __init__(self, spans=[], yspan=[], label='', style='orange', alpha=None, pos='above', size=8, center=True):
		Vspans.__init__(self, spans=spans, style=style, alpha=alpha)
		self.label = label
		self.yspan = yspan
		self.pos = pos
		self.size = size
		self.center = center

	def draw(self, plot):
		''' called when rendering plots '''
		self.set_color()

		diffspans = [[span[0], span[1]-span[0]] for span in self.spans]

		plot.ax.broken_barh(diffspans, self.yspan, facecolor=self.style, edgecolor=(1,1,1,0.5), zorder=2.0)

		if self.center: xtext = (self.spans[0][0] + self.spans[-1][-1])/2
		else: xtext = self.spans[0][0]
		if self.pos == 'above':
			xytext = [xtext, sum(self.yspan)+0.01]
		else:
			xytext = [xtext, self.yspan[0]+0.01]
		obj = plot.ax.annotate(self.label, xy=xytext, size=self.size)
		if self.center: obj.set_horizontalalignment('center')
		else: obj.set_horizontalalignment('left')

class Wedge(Entity):
	''' a vertical black bar and a little triangle pointing into the interval '''
	def __init__(self, x, scale=1, y=None, style=None):
		''' constructor
		x: left edge
		scale: scale of marker. Negative values define left-pointing markers (default:1)
		y: y position
		style: color specifier
		'''
		self.x, self.y = x, y
		self.style = 'black' if style is None else style
		self.scale = scale

	def draw(self, plot):
		''' called when rendering plots '''
		if self.y > 0: ymin, ymax = 0.5, 1
		elif self.y < 0: ymin, ymax = 0, 0.5
		else: ymin, ymax = 0, 1
		x = max(plot.xlim[0], min(plot.xlim[1], self.x))
		plot.ax.axvline(x=x, color=self.style, ymin=ymin, ymax=ymax)

		xlim = plot.ax.get_xlim()
		if (xlim[1] - xlim[0]) > 600: k = 1
		else: k = (xlim[1] - xlim[0]) / MISTAKEFUDGE

		if self.scale:
			x = x + abs(self.scale)**.5 * (self.scale)/abs(self.scale) * k
			if self.y == 0: y = 2
			else: y = self.y

			size = abs(self.scale)

			if self.scale < 0: marker = '<'
			if self.scale > 0: marker = '>'

			plot.ax.scatter([self.x], [self.y], marker=marker, color=self.style, s=25*abs(self.scale))

class Wall(Vspans):
	''' vertical black bars with little triangles pointing into the interval '''
	def __init__(self, spans=[], y=None, ylim=[0,1], style='black', wedge=1, single=False):
		''' constructor
		spans: iterable of iterables defining intervals, e.g. [[20, 45], [60, 80]] (default:[])
		y: position of triangle (default:None)
		ylim: where the black bars are relative to the plot (default:[0, 1])
		style: color specifier (default:'black')
		wedge: scale of triangle. Negative values result in left-pointing triangles (default:1)
		single: define single walls instead of pairs of walls. (default:False)
		'''
		Vspans.__init__(self, spans=spans, style=style)
		self.y = y
		self.ylim = ylim
		self.wedge = wedge
		self.single = single

	def get_y(self):
		''' get y '''
		if type(self.y) is None: return 2
		else: return self.y

	def get_ylim(self):
		''' get ylim '''
		if self.ylim is None: return [0, 1]
		elif type(self.ylim) is int or type(self.ylim) is float:
			if self.ylim > 0: return [0.5, 1]
			elif self.ylim == 0: return [0, 1]
			else: return [0, 0.5]
		else: return self.ylim

	def draw(self, plot):
		''' called when rendering plots '''
		ylim = self.get_ylim()
		for span in self.spans:
			for i, x in enumerate(span): 
				plot.ax.axvline(x=x, color='black', ymin=ylim[0], ymax=ylim[1])
				if self.wedge:
					if self.single:
						if self.wedge >= 0: marker = '>'
						else: marker = '<'
					else:
						n = i
						if self.wedge < 0: n += 1
						if n % 2: marker = '<'
						else: marker = '>'

					xlim = plot.ax.get_xlim()
					#wx = x# + (abs(4*self.wedge)**0.5 * np.sign(0.5 - i % 2))
					#FIXME: make the arrows aligned at any scale
					xlim = plot.ax.get_xlim()
					if (xlim[1] - xlim[0]) > 600: k = 1
					else: k = (xlim[1] - xlim[0]) / MISTAKEFUDGE

					wx = x - (2*(abs(self.wedge) * np.sign((i % 2) - 0.5)) - self.wedge*.5) * k

					if self.y is None: wy = 2
					else: wy = self.y

					plot.ax.scatter([wx], [wy], marker=marker, color='black', s=25*abs(self.wedge))

class HMMTOP(Vspans):
	''' Vspans but based on HMMTOP predictions '''
	def __init__(self, gseq, style='orange', alpha=None, nohmmtop=False, load=None):
		''' constructor
		gseq: sequence, if any
		style: color specifier (default:'orange')
		alpha: transparency (default:None)
		nohmmtop: don't run HMMTOP (default:False)
		load: load an existing prediction
		'''
		Vspans.__init__(self, style=style, alpha=alpha)
		self.spans = []
		self.topout = ''

		gseq = gseq.upper()
		fasta = gseq
		#no FASTA? no problem

		if load is not None:
			with open(load) as f:
				indices = self.parse_hmmtop(f.read().decode('utf-8'))
		elif nohmmtop: indices = []
		else: 
			if not fasta.startswith('>'): fasta = '>seq\n' + fasta
			p = subprocess.Popen(['hmmtop', '-if=--', '-is=pseudo', '-sf=FAS', '-pi=spred'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate(input=fasta)
			print(err.strip(), file=sys.stderr)

			indices = self.parse_hmmtop(out)

		if not indices: return

		gseq = re.sub('\*', '-', gseq)
		gseq = gseq.upper()
		if gseq.startswith('>'): seq = gseq[gseq.find('\n')+1:]
		else: seq = gseq

		seq = re.sub('[^A-Z\-]', '', seq)

		for i, c in enumerate(seq):
			if c in 'ACDEFGHIKLMNPQRSTVWY': pass
			else: 
				for j, n in enumerate(indices):
					if (n - 1) >= i: indices[j] += 1
		self.spans = [[indices[i], indices[i+1]] for i in range(0, len(indices), 2)]

	def parse_hmmtop(self, topout):
		''' parse topout into spans '''
		self.topout = topout
		if not topout: return []
		indices = re.findall('(?: IN| OUT)((?:\s*(?:[0-9]+))+)', topout.strip())[0].strip().split()
		indices = [int(i) for i in indices[1:]]
		return indices

	def add(self, *spans):
		''' add TMSs '''
		for span in spans: self.spans.append(span)

	def replace(self, *spans):
		''' replace TMSs '''
		self.delete(*spans)
		for newspan in spans: self.spans.append(newspan)

	def delete(self, *spans):
		''' delete TMSs '''
		removeme = []
		for j, oldspan in enumerate(self.spans):
			for i, newspan in enumerate(spans):
				if overlap(newspan, oldspan): 
					removeme.append(j)
					break
		for j in removeme[::-1]: 
			self.spans.pop(j)

	def extend(self, *spans):
		''' extend TMSs '''
		fuseme = {}
		appendme = []
		for i, newspan in enumerate(spans):
			notfound = 1
			for j, oldspan in enumerate(self.spans):
				if overlap(newspan, oldspan): 
					try: fuseme[i].append(j)
					except KeyError: fuseme[i] = [j]
					notfound = 0
			if notfound: appendme.append(newspan)
		#schedule fusions
		for i in fuseme:
			newtms = spans[i]
			for j in fuseme[i]: 
				newtms = union(newtms, self.spans[j])
			appendme.append(newtms)
		#remove fusion participants
		popme = set()
		for i in fuseme:
			for j in fuseme[i]:
				popme.add(j)
		for i in sorted(popme)[::-1]: self.spans.pop(i)
		#add new TMSs
		for span in appendme: 
			#if type(span) is int: self.spans.append(spans[span])
			self.spans.append(span)

class Point(Entity):
	''' a thing with an x and a y position '''
	def __init__(self, pos, marker='.', style='r', size=25):
		''' constructor 
		pos: iterable with x and y
		marker: marker at point (default:'.')
		style: color (default:'r')
		size: scale of marker (default:25)
		'''
		self.x, self.y = pos
		self.marker = marker
		self.style = style
		self.size = size

	def draw(self, plot):
		''' called when rendering plots '''
		plot.ax.scatter([self.x], [self.y], marker=self.marker, color=self.style, s=self.size)
		
class Topoprob(Curve):
	''' a curve where y is equal to the log probability that HMMTOP's predictions were accurate '''
	def __init__(self, gseq, topout, style='y', offset=0, window=19, invert_orientation=False, smooth_loops=0):
		''' constructor
		gseq: sequence
		topout: HMMTOP or HMMTOP-like output
		style: color specifier (default:'y')
		offset: shift curve this far (default:0)
		window: window size to compute entropy over (default:19)
		invert_orientation: use opposite probabilities for O/o/I/i
		smooth_loops: 0: no (default); 1: assume all O are o and I are i; 2: assume all o are O and i are I
		'''
		Curve.__init__(self, style=style)
		self.window = window
		if not topout: 
			self.X = []
			self.Y = []
			return

		if gseq.startswith('>'):
			gseq = gseq[gseq.find('\n')+1:]
			gseq = re.sub('[^A-Z\-]', '', gseq)
		gseq = gseq.upper()
		gseq = re.sub('\*', '-', gseq)
		seq = re.sub('[X\-]', '', gseq)
		seq = re.sub('[^A-Z]', '', seq)
		prev = 0

		window = 5
		aa2index = dict(zip('ACDEFGHIKLMNPQRSTVWY', range(20)))

		if invert_orientation: state2index = dict(zip('IiOoH', range(5)))
		else: state2index = dict(zip('OoIiH', range(5)))

		probs = []

		if 'HMMTOP_PSV' in os.environ and os.path.isfile(os.environ['HMMTOP_PSV']): fn = os.environ['HMMTOP_PSV']
		elif 'hmmtop.psv' in os.listdir('.'): fn = './hmmtop.psv'
		else: raise IOError('Could not find hmmtop.psv')

		skip = 1
		states = 5
		freqs = []
		with open(fn) as f:
			for l in f: 
				if skip: 
					skip -= 1
					continue
				if states:
					freqs.append([])
					for i in range(0, len(l), 8): 
						try: 
							n = float(l[i:i+8])
							freqs[-1].append(n)
						except ValueError: pass
					states -= 1
		freqmatrix = np.array(freqs)
		for i in range(freqmatrix.shape[1]): freqmatrix[:,i] /= sum(freqmatrix[:,i])
		freqmatrix = np.log2(freqmatrix/0.2)

		tags = {}
		for resi, resn in enumerate(seq): tags[resi+1] = 'X'

		essentials = re.findall('(?:OUT|IN)(?:(?:\s+[0-9]+)+)\s*$', topout)[0].strip().split()
		if essentials[0] == 'OUT': tags[0] = 'O'
		else: tags[0] = 'I'

		spans = []
		for i in range(int(essentials[1])):
			span = int(essentials[2*i + 2]), int(essentials[2*i + 2 + 1])
			spans.append(span)
			for i in range(span[0], span[1]+1):
				tags[i] = 'H'

		flipped = False
		state = 0 if tags[0] == 'O' else 'I'
		for i in sorted(tags)[1:]:
			if tags[i] == 'X':
				if state == 0: tags[i] = 'O'
				else: tags[i] = 'I'
				flipped = False
			elif tags[i] == 'H' and not flipped:
				state = not state
				flipped = True

		if smooth_loops == 1:
			for i in sorted(tags):
				if tags[i] != 'H': tags[i] = tags[i].lower()
		if smooth_loops == 2:
			for i in sorted(tags):
				if tags[i] != 'H': tags[i] = tags[i].upper()
		else:
			for span in spans:
				for n in range(1, 15+1):
					i = span[0] - n
					if (2 <= i <= len(tags)) and (tags[i] != 'H'): 
						tags[i] = tags[i].lower()
					i = span[1] + n
					if (2 <= i <= len(tags)) and (tags[i] != 'H'): 
						tags[i] = tags[i].lower()

		#s = ''
		#for resi in sorted(tags): s += '{}'.format(tags[resi])
		#print(s)

		rawfreqs = []
		for resi in sorted(tags): 
			resn = seq[resi-1]
			aa_i = aa2index[resn]
			state = tags[resi]
			s_i = state2index[state]
			#print(resi, resn, state, freqmatrix[s_i,aa_i])
			#Y.append(np.log2(freqmatrix[s_i,aa_i]/0.2))
			rawfreqs.append(freqmatrix[s_i,aa_i])
		avgfreqs = []
		for i in range(0, len(rawfreqs) - window):
			sample = rawfreqs[i:i+window]
			avgfreqs.append(sum(sample)/len(sample))
		#print(min(Y), max(Y), sum(Y)/len(Y), np.std(Y))

		if len(gseq) == len(seq):
			self.Y = np.array(avgfreqs)
			self.X = np.arange(0, len(self.Y))+window//2+offset+1
		else:
			replace = re.finditer('(-+|X+)', gseq)
			inserts = {}
			for i in replace: inserts.setdefault(i.start(), i.end()-i.start())
			first = False
			newprobs = []

			for x, p in enumerate(avgfreqs):
				if x in inserts and not first:
					first = True
					newprobs += [np.nan for y in range(inserts[x])]
					newcount = x + inserts[x]

				elif not first: newprobs.append(p)

				elif first and newcount in inserts:
					newprobs += [np.nan for y in range(inserts[newcount])]
					newcount += inserts[newcount]
				else:
					newprobs.append(p)
					newcount += 1

			self.Y = np.array(newprobs)
			self.X = np.arange(offset+1, len(self.Y)+1)+window//2

class Entropy(Curve):
	''' a curve where y is equal to the average informational content in a window about each residue '''
	def __init__(self, gseq, style='r', offset=0, window=19):
		''' constructor
		gseq: sequence
		style: color specifier (default:'r')
		offset: shift curve this far (default:0)
		window: window size to compute entropy over (default:19)
		'''
		Curve.__init__(self, style=style)
		self.window = window

		if gseq.startswith('>'):
			gseq = gseq[gseq.find('\n')+1:]
			gseq = re.sub('[^A-Z\-]', '', gseq)
		gseq = gseq.upper()
		gseq = re.sub('\*', '-', gseq)
		seq = re.sub('[X\-]', '', gseq)
		seq = re.sub('[^A-Z]', '', seq)
		prev = 0

		entropies = []
		for i in range(0, len(seq)-window):
			probs = {}
			for j in range(i, i+window):
				try: probs[seq[j]] += 1./window
				except KeyError: probs[seq[j]] = 1./window

			s = 0
			for r in probs: s += (-probs[r]*np.log2(probs[r]))
			entropies.append(s)

		self.Y = np.array(entropies)
		self.X = np.arange(0, len(self.Y))+window//2+offset+1

	def draw(self, plot):
		''' called when rendering plots '''
		self.set_color()
		plot.xlim[0] = min(plot.xlim[0], self.X[0] - self.window//2)
		plot.xlim[1] = max(plot.xlim[1], self.X[-1] + self.window//2)
		plot.ylim[0] = 2
		plot.ylim[1] = np.log2(20)

		plot.axeslabels = ['Residue number', 'Entropy (bits)']
		plot.ax.plot(self.X, self.Y, color=self.style, linewidth=1)

class Psipred(Curve):
	''' a curve where y is equal to the helix probability (red), sheet probability (yellow), and coil probability (green) as given by PSIPRED '''
	def __init__(self, pred, style=None, offset=0, window=1):
		''' constructor
		pred: PSIPREd prediction as a str
		style: color specifier (default:None)
		offset: shift x-wise (default:0)
		window: compute moving averages over this many residues (default:1)
		'''
		Curve.__init__(self, style=style)
		self.window = window

		self.parse(pred)

	def parse(self, pred):
		''' parse predictions '''
		x = []
		ycoil = []
		yhelx = []
		ystrn = []
		for l in pred.split('\n'):
			if not l.strip(): continue
			elif l.lstrip().startswith('#'): continue
			elif l.lstrip().startswith('>'): continue
			sl = l.strip().split()
			ycoil.append(float(sl[3]))
			yhelx.append(float(sl[4]))
			ystrn.append(float(sl[5]))
			x.append(int(sl[0]))

		if self.window == 1:
			self.X = np.array(x)
			self.Ycoil = np.array(ycoil)
			self.Yhelx = np.array(yhelx)
			self.Ystrn = np.array(ystrn)
		else:
			self.X = np.array(x[self.window//2:len(x)-self.window+self.window//2])
			self.Ycoil = np.array([np.mean(ycoil[i:i+self.window]) for i, span in enumerate(ycoil[:-self.window])])
			self.Yhelx = np.array([np.mean(yhelx[i:i+self.window]) for i, span in enumerate(yhelx[:-self.window])])
			self.Ystrn = np.array([np.mean(ystrn[i:i+self.window]) for i, span in enumerate(ystrn[:-self.window])])

	def draw(self, plot):
		''' called when rendering plots '''
		plot.xlim[0] = min(plot.xlim[0], self.X[0] - self.window//2)
		plot.xlim[1] = max(plot.xlim[1], self.X[-1] + self.window//2)
		plot.ylim[0] = 0
		plot.ylim[1] = 1.0
		plot.yticks = 0.2

		plot.axeslabels = ['Residue number', 'Confidence']
		plot.ax.plot(self.X, self.Ycoil, color='g', linewidth=2)
		plot.ax.plot(self.X, self.Yhelx, color='r', linewidth=2)
		plot.ax.plot(self.X, self.Ystrn, color='y', linewidth=2)

class Hydropathy(Curve):
	''' a curve where y is the moving average of the hydropathy '''
	def __init__(self, gseq, style='r', offset=0, window=19, index=HYDROPATHY):
		''' constructor
		gseq: sequence
		style: color specifier (default:'r')
		offset: shift x-wise by this amount (default:0)
		window: window size (default:19)
		'''
		Curve.__init__(self, style=style)
		self.window = window
		#Copied with barely any modification from gblast3.py

 		#sanitize the sequence, which may be a FASTA
		if gseq.startswith('>'):
			gseq = gseq[gseq.find('\n')+1:]
			gseq = re.sub('[^A-Z\-]', '', gseq)
		gseq = gseq.upper()
		gseq = re.sub('\*', '-', gseq)
		self.gseq = gseq
		seq = re.sub('[X\-]', '', gseq)
		seq = re.sub('[^A-Z]', '', seq)
		self.seq = seq
		prev = 0

		midpt = (window+1)//2
		length = len(seq)
		hydro = []
		for i in range(length-window+1):
			total = 0
			for j in range(window): total += index[seq[i+j]]
			total /= window
			hydro.append(total)

		if len(seq) == len(gseq): 
			self.Y = np.array(hydro)
			self.X = np.arange(0, len(self.Y))+window//2+offset+1

		else:
			replace = re.finditer('(-+|X+)', gseq)
			inserts = {}
			for i in replace: inserts.setdefault(i.start(), i.end()-i.start())
			first = False
			newhydro = []

			for x, h in enumerate(hydro):
				if x in inserts and not first:
					first = True
					newhydro += [np.nan for y in range(inserts[x])]
					newcount = x + inserts[x]

				elif not first: newhydro.append(h)

				elif first and newcount in inserts:
					newhydro += [np.nan for y in range(inserts[newcount])]
					newcount += inserts[newcount]
				else:
					newhydro.append(h)
					newcount += 1

			self.Y = np.array(newhydro)
			self.X = np.arange(offset+1, len(self.Y)+1)+window//2

	def draw(self, plot):
		''' called when rendering plots '''
		self.set_color()
		plot.xlim[0] = min(plot.xlim[0], self.X[0] - self.window//2)
		plot.xlim[1] = max(plot.xlim[1], self.X[-1] + self.window//2)

		plot.axeslabels = ['Residue number', 'Hydropathy (kcal/mol)']
		plot.ax.plot(self.X, self.Y, color=self.style, linewidth=1)


class What(Entity):
	''' a container with both a curve and HMMTOP vspans '''
	def __init__(self, seq, style=0, tmscolor=None, linecolor=None, color=None, nohmmtop=False, mode='hydropathy', window=None, predfile=None, topo_prob=False):
		''' constructor
		seq: sequence or input data
		style: color specifier (default:0)
		tmscolor: vspan color specifier (default:None)
		linecolor: curve color specifier (default:None)
		color: color specifier for everything (default:None)
		nohmmtop: do not run HMMTOP (default:False)
		mode: what kind of plot to generate ([hydropathy], entropy, psipred)
		window: window size for computing moving averages (default:None)
		predfile: prediction file to load for PSIPRED (default:None)
		'''
		Entity.__init__(self)
		if mode == 'psipred': self.seq = self.get_psipred_seq(seq)
		else: self.seq = seq

		self.style = style

		def choose(*options):
			for x in options[::-1]:
				if x is not None: return x
			return None

		self.tmscolor = choose(color, tmscolor)
		self.linecolor = choose(color, linecolor)

		self.entities = []
		if mode == 'entropy':
			if window is None: self.window = 19
			else: self.window = window
			self.entities.append(Entropy(seq, style=self.get_curve_color(), window=self.window))
		elif mode == 'psipred':
			if window is None: self.window = 1
			else: self.window = window
			self.entities.append(Psipred(seq, style=self.get_curve_color(), window=self.window))
		else: 
			if window is None: self.window = 19
			else: self.window = window
			self.entities.append(Hydropathy(seq, style=self.get_curve_color(), window=self.window))
		self.entities.append(HMMTOP(self.seq, style=self.get_tms_color(), nohmmtop=nohmmtop))
		if topo_prob: self.entities.append(Topoprob(seq, topout=self.entities[-1].topout, window=self.window))

	def get_psipred_seq(self, seq):
		''' obtain sequence given a prediction file '''
		s = ''
		for l in seq.split('\n'): 
			if not l.strip(): continue
			elif l.lstrip().startswith('#'): continue
			elif l.lstrip().startswith('>'): s += '\n'+ l.strip() + '\n'
			else: s += l.strip().split()[1]
		return s.strip()

	def get_title(self, showlength=True):
		''' generate a title given attributes '''
		if self.seq.startswith('>'): 
			s = self.seq[1:self.seq.find('\n')]
			if showlength: 
				truelength = len(re.sub('[^A-Z]', '', self.seq[self.seq.find('\n')+1:]))
				s += ' ({}aa)'.format(truelength)
			return s
		else: 
			PREVIEW = [8, 3]
			if len(self.seq) <= sum(PREVIEW): return '{} ({}aa)'.format(self.seq, len(self.seq))
			else: return '{}...{} ({}aa)'.format(self.seq[:8], self.seq[-3:], len(self.seq))

	def get_curve_color(self):
		''' generate a curve color '''
		if self.linecolor is None:
			if (self.style % 3) == 0: return 'red'
			elif (self.style % 3) == 1: return 'blue'
			elif (self.style % 3) == 2: return 'green'
		#elif self.linecolor.startswith('#') or self.linecolor.startswith('0x'):
		#	return hex2tuple(self.linecolor)
		else: return self.linecolor 

	def get_tms_color(self):
		''' generate a vspan color '''
		if self.tmscolor is None:
			if (self.style % 3) == 0: return 'orange'
			elif (self.style % 3) == 1: return 'cyan'
			elif (self.style % 3) == 2: return 'purple'
		#elif self.tmscolor.startswith('#') or self.tmscolor.startswith('0x'):
		#	return hex2tuple(self.tmscolor)
		else: return self.tmscolor

	def draw(self, plot):
		''' called when rendering plots '''
		#for e in self.entities: e.draw(notitle=True)
		for e in self.entities: e.draw(plot)

		#XXX maybe this goes better under Hydropathy or something
		plot.titles.append(self.get_title())

	def __iter__(self): return iter(self.entities)

def split_fasta(fastastr):
	''' split multirecord FASTAs '''
	sequences = []
	for l in fastastr.split('\n'):
		if l.startswith('>'): sequences.append(l)
		elif not l.strip(): continue
		else: 
			try: sequences[-1] += '\n' + l.strip()
			except IndexError: sequences.append('>untitled sequence\n{}'.format(l))
	return sequences

def parse_manual_tms(top, s):
	''' under construction

	will allow TMS specification using a nice mini-language '''
	tmss = []

	indices = []
	command = 'add'
	color = 'orange'
	def isint(token):
		try:
			int(token)
			return True
		except ValueError: return False
	def isverb(token):
		if token in []: pass
	for token in s.split():
		if isverb(token): pass
		if isint(token): indices.append(int(token))
	return indices

def find(entitylist, target, index=None):
	''' look for a specific kind of entity in a list of entities

	entitylist: list of entities
	target: type of entity
	index: does something (default:None)
	'''
	if index is None: anymode = 1
	else:
		anymode = 0
		try: iter(index)
		except TypeError: index = [index]
	i = 0
	out = []
	for e in entitylist:
		if type(e) is target:
			if anymode: out.append(e)
			elif i in index: out.append(e)
			i += 1
	return out

#def main(infiles, mode='hydropathy', walls=None, ):
#def what(sequences, labels=None, imgfmt='png', directory=None, filename=None, title=False, dpi=80, hide=True, viewer=None, bars=[], color='auto', offset=0, statistics=False, overwrite=False, manual_tms=None, wedges=[], ywedge=2, legend=False, window=19, silent=False, axisfont=None, tickfont=None, xticks=None, mode='hydropathy', width=None, height=None):
def what(*args, **kwargs): main(*args, **kwargs)
def main(infiles, **kwargs):
	''' generate a QUOD figure from start to finish
	will attempt to generate multi-subplot figures if multiple modes are given
	infiles: iterable of filenames

	keyword arguments:
		add_marker: add a point marker at these points
		add_region: add a region marker at these spans
		add_tms: add a TMS marker at these spans
		axis_font: font size for axis labels
		bars: add vertical bars at these x-values
		color: color everything this
		delete_tms: delete TMSs intersecting these spans
		entities: preload in this iterable of Entities
		entropy: plot entropy too
		extend_tms: extend TMSs intersecting these spans
		force_seq: force inputs to be interpreted as sequences
		load_tms: load TMS definitions from a file
		manual_tms: not implemented
		no_tms: don't run HMMTOP for these sequence IDs
		outdir: where to output figure images
		outfile: what to name figure images
		psipred: plot PSIPRED predictions too
		replace_tms: replace TMSs intersecting these spans
		tick_font: font size for tick labels
		title: manually set title
		wall: vertical black bars with wedges
		walls: vertical black bars with inward-pointing wedges
		window: window width
		x_offset: shift x-wise by this amount
		xticks: x-tick interval
	'''

	#grep -o 'kwargs\[.*:' quod.py | sed "s/kwargs\[\'//g;s/']//g;s/is not None//g;s/mode == '//g;s/'//g;s/ //g" | sort | uniq | pbcopy
	plot = Plot()

	if 'width' in kwargs: width = kwargs['width']
	else: width = None
	if 'height' in kwargs: height = kwargs['height']
	else: height = None

	plot.width = width
	plot.height = height

	if 'color' in kwargs and kwargs['color'] is not None: color = kwargs['color']
	else: color = 'auto'

	no_tms = []
	loadme = {}
	if 'load_tms' in kwargs and kwargs['load_tms']:
		n_id = 0
		n_loads = 0
		for token in kwargs['load_tms']:
			if isid(token): 
				no_tms.append(parse_id(token))
				n_id += 1
			else: 
				if n_id:
					if (no_tms[-1] % 3) == 0: style = 'orange'
					elif (no_tms[-1] % 3) == 1: style = 'cyan'
					elif (no_tms[-1] % 3) == 2: style = 'purple'
					else: style = 'orange'
				else:
					if (n_loads % 3) == 0: style = 'orange'
					elif (n_loads % 3) == 1: style = 'cyan'
					elif (n_loads % 3) == 2: style = 'purple'
					else: style = 'orange'
				x = HMMTOP(gseq='', nohmmtop=True, load=token, style=style)
				loadme[n_id] = x
				plot.add(x)
				n_loads += 1
		if n_id == 0: no_tms.append(0)

	#if 'entropy' in kwargs and kwargs['entropy']: entropy = True
	#else: entropy = False
	if 'mode' in kwargs:
		if kwargs['mode'] == 'entropy': mode = 'entropy'
		elif kwargs['mode'] == 'psipred': mode = 'psipred'
		else: mode = 'hydropathy'
	else: mode = 'hydropathy'

	if 'window' in kwargs and kwargs['window'] is not None: window = kwargs['window']
	else: window = 19

	if 'topo_prob' in kwargs and kwargs['topo_prob']: topo_prob = True
	else: topo_prob = False

	if 'no_tms' in kwargs and kwargs['no_tms']:
		for token in kwargs['no_tms']:
			if not isid(token): raise ValueError('--no-tms: Not an id: "{}"'.format(token))
			#for e in find(entities, What, parse_id(token)):
			#	for ee in find(e.entities, HMMTOP): ee.spans = []
			no_tms.append(parse_id(token))
				
	n = 0
	if 'force_seq' in kwargs and kwargs['force_seq']: 
		for seq in infiles: 
			if n in no_tms: nohmmtop = 1
			else: nohmmtop = 0
			whatkwargs = {'nohmmtop':nohmmtop, 'window':window, 'topo_prob':topo_prob}
			if color == 'auto': whatkwargs['style'] = n
			else: whatkwargs['color'] = color
			if mode == 'entropy': whatkwargs['mode'] = 'entropy'
			elif mode == 'psipred': whatkwargs['mode'] = 'psipred'
			plot.add(What(seq, **whatkwargs))

			if n in loadme: plot[-1].add(loadme[n])
			n += 1
	else:
		for fn in infiles:
			with open(fn) as f:
				for seq in split_fasta(f.read().decode('utf-8')):
					if n in no_tms: nohmmtop = 1
					else: nohmmtop = 0
					whatkwargs = {'nohmmtop':nohmmtop, 'window':window, 'topo_prob':topo_prob}
					if color == 'auto': whatkwargs['style'] = n
					else: whatkwargs['color'] = color
					if mode == 'entropy': whatkwargs['mode'] = 'entropy'
					elif mode == 'psipred': whatkwargs['mode'] = 'psipred'
					plot.add(What(seq, **whatkwargs))

					if n in loadme: plot[-1].add(loadme[n])
					n += 1

	if 'add_tms' in kwargs and kwargs['add_tms']:
		for tms in kwargs['add_tms']:
			stms = tms.split(':')
			if len(stms) >= 2: color = stms[1]
			else: color = 'orange'
			for span in stms[0].split(','):
				if len(span.split('-')) == 1: indices = [int(span)]*2
				else: indices = [int(x) for x in span.split('-')]
				spans = [indices[i:i+2] for i in range(0, len(indices), 2)]
				plot.add(HMMTOP('', style=color, alpha=None, nohmmtop=True))
				plot[-1].spans = spans

	if 'mark_residue' in kwargs and kwargs['mark_residue']:
		for marker in kwargs['mark_residue']:
			pos = []
			color = None
			index = 0
			size = 25

			motifs = []
			for i, token in enumerate(marker.split(':')):
				if isid(token) and not i: index = parse_id(token)
				elif not motifs: motifs += token.split(',')
				else: color = token

			for e in find(plot, What, index):
				for ee in find(e, Hydropathy):
					if color is None: color = ee.style
					pos = []
					for m in motifs:
						p = m.replace('X', '.').replace('-', '.?')
						for hit in re.finditer(p, ee.gseq): 
							#pos.append(hit.start())
							pos.append((hit.start() + hit.end())//2)
					done = []
					for pair in zip(ee.X, ee.Y):
						for x in pos:
							if x == pair[0]:
								if np.isnan(pair[1]): plot.add(Point([x,0], marker='o', style=color, size=size))
								else: plot.add(Point(pair, marker='o', style=color, size=size))
								done.append(x)
					for x in pos:
						if x not in done: plot.add(Point([x,0], marker='o', style=color, size=size))

	if 'add_marker' in kwargs and kwargs['add_marker']:
		for marker in kwargs['add_marker']:
			pos = []
			color = None
			index = 0
			size = 25

			for i, token in enumerate(marker.split(':')): 
				if isid(token) and not i: index = parse_id(token)
				elif isint(token[0]) and not pos:
					pos = [int(n) for n in token.split(',')]
				elif isint(token[0]) and pos:
					size = int(token)
				else: color = token

			done = []
			for e in find(plot, What, index):
				for ee in find(e, Hydropathy):
					if color is None: color = ee.style

					for pair in zip(ee.X, ee.Y):
						for x in pos:
							if x == pair[0]:
								if np.isnan(pair[1]): plot.add(Point([x,0], marker='o', style=color, size=size))
								else: plot.add(Point(pair, marker='o', style=color, size=size))
								done.append(x)
					for x in pos:
						if x not in done:
							plot.add(Point([x,0], marker='o', style=color, size=size))

	if 'extend_tms' in kwargs and kwargs['extend_tms']:
		for tms in kwargs['extend_tms']:
			stms = tms.split(':')
			index = 0
			spans = []
			color = None

			for i, token in enumerate(stms):
				if i == 0 and isid(token):
					index = parse_id(token)
				elif not spans and isspans(token): spans = parse_spans(token)
				else: color = token

			for e in find(plot, What, index):
				for ee in e.entities:
					if type(ee) is HMMTOP:
						if color is None or color == ee.style: 
							color = ee.style
							for span in spans: ee.extend(span)
						else:

							#search
							for i, oldspan in enumerate(ee.spans):
								for j, newspan in enumerate(spans):
									if overlap(oldspan, newspan): ee.delete(newspan)
							#integration
							plot.add(HMMTOP('', style=color))
							plot[-1].spans = spans

	#if an existing TMS on sequence +x overlaps with a TMS defined by delete_tms, erase it
	if 'delete_tms' in kwargs and kwargs['delete_tms']:
		for tms in kwargs['delete_tms']:
			stms = tms.split(':')

			index = 0
			spans = []

			for i, token in enumerate(stms):
				if isid(token) and not i: index = parse_id(token)
				elif not spans and isspans(token): spans = parse_spans(token)
				else: color = token

			for e in find(plot, What, index):
				for ee in find(e.entities, HMMTOP): ee.delete(*spans)

	if 'replace_tms' in kwargs and kwargs['replace_tms']:
		for tms in kwargs['replace_tms']:
			stms = tms.split(':')
			index = 0
			spans = []
			color = None

			for i, token in enumerate(stms):
				if isid(token) and not i: index = parse_id(token)
				elif not spans and isspans(token): spans = parse_spans(token)
				else: color = token

			for e in find(plot, What, index):
				for ee in find(e.entities, HMMTOP):
					if color is None or color == ee.style: ee.replace(*spans)
					else:
						ee.delete(*spans)
						plot.add(HMMTOP('', style=color))
						plot[-1].spans = spans

	if 'add_region' in kwargs and kwargs['add_region']:
		for region in kwargs['add_region']:
			tokens = ['']
			skip = 0
			nosplit = 0
			for i, c in enumerate(region):
				if skip:
					skip -= 1
					continue
				if c == '\\': 
					try:
						if region[i+1] == ':': 
							nosplit = 1
							tokens[-1] += c
					except IndexError: tokens[-1].append(c)
				elif nosplit:
					tokens[-1] += c
					nosplit -= 1
				elif c == ':': tokens.append('')
				else: tokens[-1] += c
			spans = []
			y = None
			label = None
			color = None
			size = None
			alpha = None

			for token in tokens:
				if isspans(token) and not spans: spans = parse_spans(token)
				elif isfloat(token) and y is None: y = float(token)
				elif label is None: label = token.replace('\\', '')
				elif color is None: color = token
				elif isint(token) and size is None: size = int(token)

			for i in range(len(spans)): spans[i][1] = spans[i][1]# - spans[i][0]

			if y is None: y = -2
			if label is None: label = 'untitled region'
			if size is None: size = 8

			plot.add(Region(spans, [y-0.15, 0.15], label, style=color, size=size))
			#for token in re.split(r':', region): print(token)

	if 'xticks' in kwargs and kwargs['xticks'] is not None: plot.xticks = kwargs['xticks']
	else: plot.xticks = None

	if 'bars' in kwargs and kwargs['bars'] is not None: 
		[plot.add(wall) for wall in parse_walls(kwargs['bars'], wedge=0)]

	if 'dpi' in kwargs: dpi = kwargs['dpi']
	else: dpi = 80

	if 'imgfmt' in kwargs: imgfmt = kwargs['imgfmt']
	else: imgfmt = 'png'

	if 'walls' in kwargs and kwargs['walls'] is not None: 
		[plot.add(wall) for wall in parse_walls(kwargs['walls'])]

	if 'wall' in kwargs and kwargs['wall'] is not None:
		[plot.add(wall) for wall in parse_walls(kwargs['wall'], single=1)]

	if 'outdir' in kwargs and kwargs['outdir']: prefix = kwargs['outdir']
	else: prefix = ''

	if 'outfile' in kwargs and kwargs['outfile']: 
		if prefix: outfile = '{}/{}'.format(prefix, kwargs['outfile'])
		else: outfile = '{}'.format(kwargs['outfile'])
	else: outfile = None

	if 'quiet' in kwargs: quiet = kwargs['quiet']
	else: quiet = True

	if 'viewer' in kwargs: viewer = kwargs['viewer']
	else: viewer = None

	if 'title' in kwargs and kwargs['title'] is not None: title = kwargs['title']
	else: title = None

	if 'x_offset' in kwargs and kwargs['x_offset'] is not None: x_offset = kwargs['x_offset']
	else: x_offset = 0
	plot.x_offset = x_offset

	if 'manual_tms' in kwargs and kwargs['manual_tms'] is not None: 
		pass #entities += parse_manual_tms(kwargs['manual_tms'])

	if 'entities' in kwargs and kwargs['entities'] is not None: 
		for e in kwargs['entities']: plot.add(e)

	if 'grid' in kwargs and kwargs['grid']: plot.grid = True

	plot.render()

	if 'axis_font' in kwargs and kwargs['axis_font'] is not None:
		plot.ax.set_xlabel(plot.ax.get_xlabel(), fontsize=kwargs['axis_font'])
		plot.ax.set_ylabel(plot.ax.get_ylabel(), fontsize=kwargs['axis_font'])

	if 'tick_font' in kwargs and kwargs['tick_font'] is not None:
		plot.ax.tick_params(labelsize=kwargs['tick_font'])

	if title is not None: plot.ax.set_title(title)
	if outfile is None: outfile = sanitize(plot.ax.get_title())

	if outfile.lower().endswith('.{}'.format(imgfmt.lower())): pass
	elif len(os.path.splitext(outfile)) == 1: outfile += '.{}'.format(imgfmt)
	elif len(os.path.splitext(outfile)[-1]) not in [3, 4]: 
		outfile += '.{}'.format(imgfmt)
	elif os.path.splitext(outfile)[-1].lower() != imgfmt.lower():
		outfile += '.{}'.format(imgfmt)

	if prefix not in outfile: outfile = '{}/{}'.format(prefix, outfile)
	plot.save(outfile, dpi=dpi, format=imgfmt)

	if not quiet:
		if viewer is None:
			if sys.platform.startswith('linux'): viewer = 'xdg-open'
			elif sys.platform.startswith('darwin'): viewer = 'open'
			else: print('[WARNING]: Unknown platform; Email khendarg@ucsd.edu with your OS and default generic file-opening program', file=sys.stderr)
		os.system('{} "{}"'.format(viewer, outfile))
		

def parse_walls(strlist, wedge=1, single=False):
	'''
	turns a list of wall specifications into Wall() objects

	strlist: a list of string with syntax "spanslike1,spanslike2,...:float:float" where spanslike is a comma-separated list of ranges, the first float defines the y position of the marker, and the second float defines the size/direction of the marker(s)
	wedge: another way to set marker size/direction
	single: produce single walls instead of double walls
	'''
	if not strlist: return None
	out = []
	for wall in strlist:
		if type(wall) is str:
			tokens = wall.split(':')
			if len(tokens) == 1: tokens.append(None)
			if len(tokens) == 2: tokens.append(None)

			if tokens[1] is None: y = None
			elif len(tokens[1]) == 0: y = None
			else: y = float(tokens[1])

			if y is None or y == 0: ylim = [0, 1]
			elif y > 0: ylim = [0.5, 1]
			else: ylim = [0, 0.5]

			wedge = wedge if tokens[2] is None else float(tokens[2])

			spans = parse_spans(tokens[0])

			out.append(Wall(spans, y=y, ylim=ylim, wedge=wedge, single=single))
		elif type(wall) is int:
			out.append(Wall([[wall, wall]], wedge=wedge, single=single))
		else: raise ValueError('Unsupported strlist format')
	return out

if __name__ == '__main__': 
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('infile', nargs='*', default=['/dev/stdin'], help='sequence files to read in')
	parser.add_argument('-a', metavar='viewer', help='Viewer to be used for opening graphs')
	parser.add_argument('-b', '--bars', nargs='+', type=int, help='Draws vertical bars at these positions')
	#-B boxes
	parser.add_argument('-c', '--color', metavar='color', default='auto', help='Colors EVERYTHING this color')
	parser.add_argument('-d', metavar='outdir', help='Directory to store graphs in (recommended only with autogenerated filenames)')
	parser.add_argument('-e', '--entropy', action='store_true', help='Plot entropy instead of hydropathy')
	#-f filename
	#-i flags file
	parser.add_argument('-l', '--title', metavar='graph_title', help='Label graph with a specific title')
	#parser.add_argument('-m', '--mode', metavar='TMSs', nargs='+', help='something something something')
	#CHECKME: is this used for anything automatic?
	parser.add_argument('-o', metavar='outfile', help='Filename of graph, relative to the argument of -d if present and as normally interpreted otherwise')
	parser.add_argument('-q', action='store_true', help='"quiet" mode, disables automatic opening of graphs')
	parser.add_argument('-r', metavar='resolution', type=int, help='Resolution of graph in dpi. The default is 80dpi, which is suitable for viewing on monitors. Draft-quality images should be about 300dpi, and publication-quality images need to be 600 or 1200dpi depending on the journal.')
	parser.add_argument('-s', action='store_true', help='Force inputs to be interpreted as sequences (this is no longer a default behavior for infile args matching /[A-Z]+/')
	parser.add_argument('-sp', action='store_true', help='Force inputs to be interpreted as PSIPRED predictions')
	parser.add_argument('-t', metavar='format', default='png', help='Format of graph (\033[1mpng\033[0m, eps, jpeg, jpg, pdf, pgf, ps, raw, rgba, svg, svgz, tif, tiff')
	parser.add_argument('-v', action='store_true', help='Verbose output. Enables warnings and generally makes things messier')
	parser.add_argument('-w', '--walls', metavar='x0-x(:scale(:y))', nargs='+', help='Draws bounds around sequences and such, e.g. "20-45,60-80:0.5:2" draws 2x-scaled markers around the interval [20, 45] and [60, 80] at y=0.5 and "36-60" draws default-scaled (1x) markers around the interval [36, 60] at the default y (y=2)')
	parser.add_argument('-W', '--wall', metavar='x(:scale(:y))', nargs='+', help='Draws a single wall for each specification, allowing left-pointing intervals with negative scale values. See -w/--walls for more information on syntax.')
	parser.add_argument('--axis-font', metavar='size', type=int, help='Axis label size (pt)')
	parser.add_argument('--grid', action='store_true', help='Show grid (default:off)')
	parser.add_argument('--height', metavar='height', type=float, help='Plot height in inches (default:5.5)')
	parser.add_argument('--mode', default='hydropathy', help='mode to run QUOD in (\033[1mhydropathy\033[0m, entropy)')
	parser.add_argument('--tick-font', metavar='size', type=int, help='Tick label size')
	parser.add_argument('--topo-prob', action='store_true', help='Plot average HMMTOP log2 emission probabilities')
	parser.add_argument('--viewer', metavar='viewer', default=None, help='Viewer to be used for opening plots')
	parser.add_argument('--width', metavar='width', type=float, help='Plot width in inches (default:dynamic)')
	parser.add_argument('--window', metavar='windowsize', type=int, default=19, help='Window size for hydropathy')
	parser.add_argument('--x-offset', metavar='init_resi', default=0, type=int, help='Sets starting x-value')
	parser.add_argument('--xticks', default=None, type=int, help='X tick spacing')

	parser.add_argument('-ar', '--add-region', metavar='x0-x1(:color)(:"label")', nargs='+')

	parser.add_argument('-am', '--add-marker', metavar='(+id):x1,x2,x3,...xn(:color)', nargs='+', help='Adds a circular marker at the specified positions on the hydropathy curve of the specified sequence')
	parser.add_argument('-mr', '--mark-residue', metavar='(+id):resn1,resn2,resn3(:color)', nargs='+', help='Adds circular markers on the specified curve on the hydropathy curve of the specified sequence')

	parser.add_argument('-at', '--add-tms', metavar='x0-x1(:color)', nargs='+', help='Adds TMSs to plot, e.g. 40-60:red 65-90:blue to add a red TMS covering residues 40 to 60 and a blue one covering residues 65 to 90 (default color: orange)')
	parser.add_argument('-dt', '--delete-tms', metavar='(+id):x0-x1,x0-x1', nargs='+', help='Deletes TMSs on plot. Use +id to specify which sequence\'s TMSs should be deleted, e.g. "+0:5-25,30-40" to delete overlapping TMSs predicted for the first sequence and "+2:60-80:red,+4:70-95:blue" to delete overlapping TMSs predicted for the third and fifth sequence. +id specification propagates rightward and defaults to +0.')
	parser.add_argument('-et', '--extend-tms', metavar='(+id):x0-x1(:color)', nargs='+', help='Extends TMSs on plot. Use +id to specify which sequence\'s TMSs should be extended, e.g. "+0:5-25,30-40" to extend TMSs predicted for the first sequence to include residues 5-25 and 30-40 without changing colors and "+2:60-80:red,+4:70-95:blue" to extend TMSs predicted for the third sequence to include residues 60 to 80 and recolor them red and TMSs predicted for the fifth sequence to include residues 70 to 95 and recolor them blue. +id specification propagates rightward and defaults to +0.')
	parser.add_argument('-lt', '--load-tms', metavar='(+id1) fn1', nargs='+', help='Loads TMS definitions from HMMTOP output stored in a file')
	parser.add_argument('-nt', '--no-tms', metavar='(+id)', nargs='+', help='Erases all TMSs for specified sequences. Applies early, so other TMS operations will override this effect.')
	parser.add_argument('-rt', '--replace-tms', metavar='(+id):x0-x1(:color)', nargs='+', help='Replaces TMSs on plot. Use +id to specify which sequence\'s TMSs should be replaced, e.g. "+0:5-25,30-40" to replace overlapping TMSs predicted for the first sequence with new TMSs spanning 5-25 and 30-40 without changing colors and "+2:60-80:red +4:70-95:blue" to replace overlapping TMSs predicted for the third sequence with a new TMS spanning residues 60 to 80 and recolor them red and overlapping TMSs predicted for the fifth sequence with a new TMS spanning residues 70 to 95 and recolor it blue. +id specification propagates rightward and defaults to +0.')

	args = parser.parse_args()

	if args.v: VERBOSITY = 1
	if not VERBOSITY:
		import warnings
		warnings.filterwarnings('ignore', '.')

	#maybe later
	#if len(args.manual_tms) == 1: tms_script = args.manual_tms[0]
	#else: 
	#	tms_script = ''
	#	for cmd in args.manual_tms:
	#		tms_script += cmd + ' '
	#tms_script = tms_script.strip()

	mode = args.mode
	if args.entropy: mode = 'entropy'
	elif args.sp: mode = 'psipred'

	main(args.infile, mode=mode, walls=args.walls, wall=args.wall, bars=args.bars, dpi=args.r, imgfmt=args.t, force_seq=args.s, outdir=args.d, outfile=args.o, color=args.color, title=args.title, quiet=args.q, viewer=args.a, axis_font=args.axis_font, width=args.width, height=args.height, x_offset=args.x_offset, add_tms=args.add_tms, delete_tms=args.delete_tms, extend_tms=args.extend_tms, replace_tms=args.replace_tms, no_tms=args.no_tms, tick_font=args.tick_font, add_marker=args.add_marker, mark_residue=args.mark_residue, add_region=args.add_region, xticks=args.xticks, load_tms=args.load_tms, entropy=args.entropy, window=args.window, grid=args.grid, topo_prob=args.topo_prob)

