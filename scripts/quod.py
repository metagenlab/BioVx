#!/usr/bin/env python2
#QUOD2: Questionable Utility of Doom 2
#most code copied from gblast3.py (which was written by Vamsee Reddy and Vasu Pranav Sai Iddamsetty) except where noted
#-Kevin Hendargo
from __future__ import division, print_function, division
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import os, subprocess, re, sys
import argparse
import numpy as np

from Bio import SeqIO

VERBOSITY = 0

def isspans(token):
	if re.match(r'^[0-9]+(?:-[0-9]+)(?:,[0-9]+(?:-[0-9]+))*$', token): return True
	else: return False

def parse_spans(token):
	spans = []
	for span in token.split(','):
		indices = span.split('-')
		if len(indices) == 0: continue
		elif len(indices) == 1: spans.append([int(indices)] * 2)
		elif len(indices) >= 2: spans.append([int(index) for index in indices[:2]])
	return spans

def isid(token):
	if re.match('^\+[0-9]+$', token): return True
	else: return False

def parse_id(token):
	if isid(token): return int(token)
	else: return None

def isint(token):
	try:
		int(token)
		return True
	except ValueError: return False

def isfloat(token):
	try:
		float(token)
		return True
	except ValueError: return False

def warn(*args): print('[WARNING]:', *args, file=sys.stderr)

def sanitize(s): return re.sub('/', '', s)

def roundto(x, n):
	return (x//5)*5

def overlap(span1, span2):
	if span1[0] <= span2[0] <= span1[-1]: return True
	elif span1[0] <= span2[1] <= span1[-1]: return True
	elif span2[0] <= span1[0] <= span2[-1]: return True
	elif span2[0] <= span1[1] <= span2[-1]: return True
	else: return False

def union(span1, span2):
	if not overlap(span1, span2): raise NotImplementedError
	else: return [min(span1[0], span2[0]), max(span1[-1], span2[-1])]

def hex2tuple(s):
	if s.startswith('#'): s = s[1:]
	elif s.startswith('0x'): s = s[2:]

	if len(s) == 3: s = [2*s[i] for i in range(len(s))]
	if len(s) == 1: s *= 6
	
	l = [int(s[i:i+2], 16)/255. for i in range(0, len(s), 2)]
	return tuple(l)

class Plot(object):
	def __init__(self): 
		self.fig = Figure()
		self.canvas = FigureCanvas(self.fig)

		self.ax = self.fig.add_subplot(111)
		self.width, self.height = None, None
		self.x_offset = 0
		self.xlim = [0, 20]
		self.ylim = [-3, 3]
		self.xticks, self.yticks = None, 1

		self.axeslabels = ['X', 'Y']
		self.titles = []

	def render(self):
		self.ax.set_xlim(left=self.xlim[0], right=self.xlim[1])
		self.ax.set_ylim(bottom=self.ylim[0], top=self.ylim[1])
		xlim = self.ax.get_xlim()
		ylim = self.ax.get_ylim()

		if self.x_offset:
			self.ax.set_xlim([x + self.x_offset for x in xlim])

		xlim = self.ax.get_xlim()
		ylim = self.ax.get_ylim()

		maxl = xlim[1] - xlim[0]
		if self.xticks is None:
			self.xticks = roundto(maxl // 5, 5)

		self.ax.set_xticks(np.arange(xlim[0], xlim[1]+1, self.xticks))
		self.ax.set_yticks(np.arange(ylim[0], ylim[1]+1, self.yticks))

		#TODO: allow centimeters or something
		#if self.width is None: self.width = (0.0265) * (maxl) if maxl > 600 else 15.
		if self.width is None: self.width = (0.0265) * (maxl) if maxl > 200 else 15.
		if self.height is None: self.height = 5.5
		self.fig.set_figwidth(self.width)
		self.fig.set_figheight(self.height)
		#self.fig.set_tight_layout(True)
		self.ax.set_xlabel(self.axeslabels[0])
		self.ax.set_ylabel(self.axeslabels[1])

		#only when plotting hydropathy! FIXME
		self.ax.axhline(y=0, color='black', linewidth=0.5)

		title = ''
		for t in self.titles: 
			if t is not None: title += '{}, '.format(t)
		self.ax.set_title(title[:-2])

	def save(self, filename, dpi=80, format='png'):
		self.fig.savefig(filename, dpi=dpi, format=format)

class Entity(object):
	def __init__(self): pass

	def draw(self, plot): raise NotImplementedError

class Curve(Entity):
	def __init__(self, X=[], Y=[], style='auto'):
		Entity.__init__(self)
		self.X = X
		self.Y = Y
		self.style = style

	def draw(self, plot): 
		if len(Y) and not len(X): 
			plot.ax.plot(Y, style)
		else: plot.ax.plot(X, Y, style)

class Vspans(Entity):
	def __init__(self, spans=[], style='orange', alpha=None):
		Entity.__init__(self)
		self.spans = spans
		self.style = style
		self.alpha = alpha

	def get_alpha(self):
		if type(self.style) is tuple and len(self.style) == 4: return self.style[3]
		elif self.alpha is not None: return self.alpha
		else: return 0.25

	def draw(self, plot):
		for span in self.spans:
			plot.ax.axvspan(span[0], span[1], facecolor=self.style, alpha=self.get_alpha())

class Region(Vspans):
	def __init__(self, spans=[], yspan=[], label='', style='orange', alpha=None, pos='above', size=8):
		Vspans.__init__(self, spans=spans, style=style, alpha=alpha)
		self.label = label
		self.yspan = yspan
		self.pos = pos
		self.size = size

	def draw(self, plot):
		plot.ax.broken_barh(self.spans, self.yspan, facecolor=self.style, edgecolor=(1,1,1,0.5), zorder=2.0)

		if self.pos == 'above':
			xytext = [self.spans[0][0], sum(self.yspan)+0.01]
		else:
			xytext = [self.spans[0][0], self.yspan[0]+0.01]
		plot.ax.annotate(self.label, xy=xytext, size=self.size)

class Wall(Vspans):
	def __init__(self, spans=[], y=None, ylim=[0,1], style='black', wedge=1, single=False):
		Vspans.__init__(self, spans=spans, style=style)
		self.y = y
		self.ylim = ylim
		self.wedge = wedge
		self.single = single

	def get_y(self):
		if type(self.y) is None: return 2
		else: return self.y

	def get_ylim(self):
		if self.ylim is None: return [0, 1]
		elif type(self.ylim) is int:
			if self.ylim > 0: return [0.5, 1]
			elif self.ylim == 0: return [0, 1]
			else: return [0, 0.5]
		else: return self.ylim

	def draw(self, plot):
		ylim = self.get_ylim()
		for span in self.spans:
			for i, x in enumerate(span): 
				plot.ax.axvline(x=x, color='black', ymin=ylim[0], ymax=ylim[1])
				if self.wedge:
					if self.single:
						if self.wedge >= 0: marker = '>'
						else: marker = '<'
					else:
						if i % 2: marker = '<'
						else: marker = '>'

					#wx = x# + (abs(4*self.wedge)**0.5 * np.sign(0.5 - i % 2))
					#FIXME: make the arrows aligned at any scale
					wx = x - 2*(abs(self.wedge) * np.sign((i % 2) - 0.5)) - self.wedge*.5

					if self.y >= 0: wy = 2
					else: wy = -2

					plot.ax.scatter([wx], [wy], marker=marker, color='black', s=25*abs(self.wedge))

class HMMTOP(Vspans):
	def __init__(self, gseq, style='orange', alpha=None):
		Vspans.__init__(self, style=style, alpha=alpha)
		self.spans = []

		fasta = gseq
		#no FASTA? no problem
		if not fasta.strip(): return

		if not fasta.startswith('>'): fasta = '>seq\n' + fasta
		p = subprocess.Popen(['hmmtop', '-if=--', '-is=pseudo', '-sf=FAS', '-pi=spred'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = p.communicate(input=fasta)
		print(err.strip(), file=sys.stderr)

		indices = re.findall('(?: IN| OUT)((?:\s*(?:[0-9]+))+)', out.strip())[0].strip().split()
		indices = [int(i) for i in indices[1:]]

		if not indices: return

		if gseq.startswith('>'): seq = gseq[gseq.find('\n')+1:]
		else: seq = gseq

		seq = re.sub('[^A-Z\-]', '', seq)

		for i, c in enumerate(seq):
			if c in 'ACDEFGHIKLMNPQRSTVWY': pass
			else: 
				for j, n in enumerate(indices):
					if (n - 1) >= i: indices[j] += 1
		self.spans = [[indices[i], indices[i+1]] for i in range(0, len(indices), 2)]

	def add(self, *spans):
		for span in spans: self.spans.append(span)

	def replace(self, *spans):
		self.delete(*spans)
		for newspan in spans: self.spans.append(newspan)

	def delete(self, *spans):
		removeme = []
		for j, oldspan in enumerate(self.spans):
			for i, newspan in enumerate(spans):
				if overlap(newspan, oldspan): 
					removeme.append(j)
					break
		for j in removeme[::-1]: 
			self.spans.pop(j)

	def extend(self, *spans):
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
	def __init__(self, pos, marker='.', style='r', size=25):
		self.x, self.y = pos
		self.marker = marker
		self.style = style
		self.size = size

	def draw(self, plot):
		plot.ax.scatter([self.x], [self.y], marker=self.marker, color=self.style, s=self.size)
		
	
class Hydropathy(Curve):

	def __init__(self, gseq, style='r', offset=0, window=19):
		Curve.__init__(self, style=style)
		self.window = window
		#Copied with barely any modification from gblast3.py

 		#sanitize the sequence, which may be a FASTA
		if gseq.startswith('>'):
			gseq = gseq[gseq.find('\n')+1:]
			gseq = re.sub('[^A-Z\-]', '', gseq)
		seq = re.sub('[X\-]', '', gseq)
		seq = re.sub('[^A-Z]', '', seq)
		prev = 0
		index = {'G':(-0.400,0.48), \
             'I':(4.500,1.38), \
             'S':(-0.800,-0.18), \
			 'Q':(-3.500,-0.85), \
             'E':(-3.500,-0.74), \
             'A':(1.800,0.62), \
             'M':(1.900,0.64), \
             'T':(-0.700,-0.05), \
             'Y':(-1.300,0.26), \
             'H':(-3.200,-0.4), \
             'V':(4.200,1.08), \
             'F':(2.800,1.19), \
             'C':(2.500,0.29), \
             'W':(-0.900,0.81), \
             'K':(-3.900,-1.5), \
             'L':(3.800,1.06), \
             'P':(-1.600,0.12), \
             'N':(-3.500,-0.78), \
             'D':(-3.500,-0.90), \
             'R':(-4.500,-2.53), \
             'U':(0,0), \
             'B':(-3.500,-0.84), \
             'J':(-3.500,-0.80), \
             'Z':(4.150,1.22) \
            }

		midpt = (window+1)//2
		length = len(seq)
		hydro = []
		for i in range(length-window+1):
			total = 0
			for j in range(window): total += index[seq[i+j]][0]
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
		plot.xlim[0] = min(plot.xlim[0], self.X[0] - self.window//2)
		plot.xlim[1] = max(plot.xlim[1], self.X[-1] + self.window//2)

		plot.axeslabels = ['Residue number', 'Hydropathy (kcal/mol)']
		plot.ax.plot(self.X, self.Y, color=self.style, linewidth=1)

class What(Entity):
	def __init__(self, seq, style=0, tmscolor=None, linecolor=None, color=None):
		Entity.__init__(self)
		self.seq = seq
		self.style = style

		def choose(*options):
			for x in options[::-1]:
				if x is not None: return x
			return None

		self.tmscolor = choose(color, tmscolor)
		self.linecolor = choose(color, linecolor)

		self.entities = []
		self.entities.append(Hydropathy(seq, style=self.get_curve_color()))
		self.entities.append(HMMTOP(seq, style=self.get_tms_color()))

	def get_title(self, showlength=True):
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
		if self.linecolor is None:
			if (self.style % 3) == 0: return 'red'
			elif (self.style % 3) == 1: return 'blue'
			elif (self.style % 3) == 2: return 'green'
		#elif self.linecolor.startswith('#') or self.linecolor.startswith('0x'):
		#	return hex2tuple(self.linecolor)
		else: return self.linecolor 

	def get_tms_color(self):
		if self.tmscolor is None:
			if (self.style % 3) == 0: return 'orange'
			elif (self.style % 3) == 1: return 'cyan'
			elif (self.style % 3) == 2: return 'purple'
		#elif self.tmscolor.startswith('#') or self.tmscolor.startswith('0x'):
		#	return hex2tuple(self.tmscolor)
		else: return self.tmscolor

	def draw(self, plot):
		#for e in self.entities: e.draw(notitle=True)
		for e in self.entities: e.draw(plot)

		#XXX maybe this goes better under Hydropathy or something
		plot.titles.append(self.get_title())

def split_fasta(fastastr):
	sequences = []
	for l in fastastr.split('\n'):
		if l.startswith('>'): sequences.append(l)
		elif not l.strip(): continue
		else: 
			try: sequences[-1] += '\n' + l.strip()
			except IndexError: sequences.append('>untitled sequence\n{}'.format(l))
	return sequences

def parse_manual_tms(top, s):
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

def find(entitylist, target, index=None):
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
def main(infiles, **kwargs):
	plot = Plot()
	entities = []

	if 'width' in kwargs: width = kwargs['width']
	if 'height' in kwargs: height = kwargs['height']

	plot.width = width
	plot.height = height

	if 'color' in kwargs and kwargs['color'] is not None: color = kwargs['color']
	else: color = 'auto'

	n = 0
	if 'force_seq' in kwargs and kwargs['force_seq']: 
		for seq in infiles: 
			if color == 'auto': entities.append(What(seq, style=n))
			else: entities.append(What(seq, color=color))
			n += 1
	else:
		for fn in infiles:
			with open(fn) as f:
				for seq in split_fasta(f.read().decode('utf-8')):
					if color == 'auto': entities.append(What(seq, style=n))
					else: entities.append(What(seq, color=color))
					n += 1

	if 'no_tms' in kwargs and kwargs['no_tms']:
		target = None
		for token in kwargs['no_tms']:
			if not isid(token): raise ValueError('--no-tms: Not an id: "{}"'.format(token))
			for e in find(entities, What, parse_id(token)):
				for ee in find(e.entities, HMMTOP): ee.spans = []
				
	if 'add_tms' in kwargs and kwargs['add_tms']:
		for tms in kwargs['add_tms']:
			stms = tms.split(':')
			if len(stms) >= 2: color = stms[1]
			else: color = 'orange'
			for span in stms[0].split(','):
				if len(span.split('-')) == 1: indices = [int(span)]*2
				else: indices = [int(x) for x in span.split('-')]
				spans = [indices[i:i+2] for i in range(0, len(indices), 2)]
				entities.append(HMMTOP('', style=color, alpha=None))
				entities[-1].spans = spans

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
			for e in find(entities, What, index):
				for ee in find(e.entities, Hydropathy):
					if color is None: color = ee.style

					for pair in zip(ee.X, ee.Y):
						for x in pos:
							if x == pair[0]:
								if np.isnan(pair[1]): entities.append(Point([x,0], marker='o', style=color, size=size))
								else: entities.append(Point(pair, marker='o', style=color, size=size))
								done.append(x)
					for x in pos:
						if x not in done:
							entities.append(Point([x,0], marker='o', style=color, size=size))

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

			for e in find(entities, What, index):
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
							entities.append(HMMTOP('', style=color))
							entities[-1].spans = spans

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

			for e in find(entities, What, index):
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

			for e in find(entities, What, index):
				for ee in find(e.entities, HMMTOP):
					if color is None or color == ee.style: ee.replace(*spans)
					else:
						ee.delete(*spans)
						entities.append(HMMTOP('', style=color))
						entities[-1].spans = spans

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

			for i in range(len(spans)): spans[i][1] = spans[i][1] - spans[i][0]

			if y is None: y = -2
			if label is None: label = 'untitled region'
			if size is None: size = 8

			entities.append(Region(spans, [y-0.15, 0.15], label, style=color, size=size))
			#for token in re.split(r':', region): print(token)

	if 'bars' in kwargs and kwargs['bars'] is not None: 
		[entities.append(wall) for wall in parse_walls(kwargs['bars'], wedge=0)]

	if 'dpi' in kwargs: dpi = kwargs['dpi']
	else: dpi = 80

	if 'imgfmt' in kwargs: imgfmt = kwargs['imgfmt']
	else: imgfmt = 'png'

	if 'walls' in kwargs and kwargs['walls'] is not None: 
		[entities.append(wall) for wall in parse_walls(kwargs['walls'])]

	if 'wall' in kwargs and kwargs['wall'] is not None:
		[entities.append(wall) for wall in parse_walls(kwargs['wall'], single=1)]

	if 'outdir' in kwargs and kwargs['outdir']: prefix = kwargs['outdir']
	else: prefix = ''

	if 'outfile' in kwargs and kwargs['outfile']: 
		if prefix: outfile = '{}/{}'.format(prefix, kwargs['outfile'])
		else: outfile = '{}'.format(kwargs['outfile'])
	else: outfile = None

	if 'quiet' in kwargs: quiet = kwargs['quiet']
	else: quiet = None

	if 'viewer' in kwargs: viewer = kwargs['viewer']
	else: viewer = None

	if 'title' in kwargs and kwargs['title'] is not None: title = kwargs['title']
	else: title = None

	if 'x_offset' in kwargs and kwargs['x_offset'] is not None: x_offset = kwargs['x_offset']
	else: x_offset = 0
	plot.x_offset = x_offset

	if 'manual_tms' in kwargs and kwargs['manual_tms'] is not None: 
		pass #entities += parse_manual_tms(kwargs['manual_tms'])

	for e in entities: e.draw(plot)

	plot.render()

	if 'axis_font' in kwargs and kwargs['axis_font'] is not None:
		plot.ax.set_xlabel(plot.ax.get_xlabel(), fontsize=kwargs['axis_font'])
		plot.ax.set_ylabel(plot.ax.get_ylabel(), fontsize=kwargs['axis_font'])

	if 'tick_font' in kwargs and kwargs['tick_font'] is not None:
		plot.ax.tick_params(labelsize=kwargs['tick_font'])

	if title is not None: plot.ax.set_title(title)
	if outfile is None: outfile = sanitize(plot.ax.get_title())

	if len(os.path.splitext(outfile)) == 1: outfile += '.{}'.format(imgfmt)
	elif len(os.path.splitext(outfile)[-1]) not in [3, 4]: 
		outfile += '.{}'.format(imgfmt)
	elif os.path.splitext(outfile)[-1].lower() != imgfmt.lower():
		outfile += '.{}'.format(imgfmt)

	plot.save(outfile, dpi=dpi, format=imgfmt)

	if not quiet:
		if viewer is None:
			if sys.platform.startswith('linux'): viewer = 'xdg-open'
			elif sys.platform.startswith('darwin'): viewer = 'open'
			else: print('[WARNING]: Unknown platform; Email khendarg@ucsd.edu with your OS and default generic file-opening program', file=sys.stderr)
		os.system('{} "{}"'.format(viewer, outfile))
		

def parse_walls(strlist, wedge=1, single=False):
	#turns a list of wall specifications into Wall() objects
	if not strlist: return None
	out = []
	for wall in strlist:
		if type(wall) is str:
			coords = [int(x) for x in wall.split(',')]
			if len(coords) == 1: coords.append(5)
			if len(coords) == 2: coords.append(None)

			span = [coords[0], coords[0]+coords[1]]
			out.append(Wall([span], y=coords[2], ylim=coords[2], wedge=wedge, single=single))
		elif type(wall) is int:
			out.append(Wall([[wall, wall]], wedge=wedge, single=single))
	return out

if __name__ == '__main__': 
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('infile', nargs='*', default=['/dev/stdin'], help='sequence files to read in')
	parser.add_argument('-a', metavar='viewer', help='Viewer to be used for opening graphs')
	parser.add_argument('-b', '--bars', nargs='+', type=int, help='Draws vertical bars at these positions')
	#-B boxes
	parser.add_argument('-c', '--color', metavar='color', default='auto', help='Colors EVERYTHING this color')
	parser.add_argument('-d', metavar='outdir', help='Directory to store graphs in (recommended only with autogenerated filenames)')
	#-e entropy
	#-f filename
	#-i flags file
	parser.add_argument('-l', '--title', metavar='graph_title', help='Label graph with a specific title')
	parser.add_argument('-m', '--manual-tms', metavar='TMSs', nargs='+', help='something something something')
	#CHECKME: is this used for anything automatic?
	parser.add_argument('-o', metavar='outfile', help='Filename of graph, relative to the argument of -d if present and as normally interpreted otherwise')
	parser.add_argument('-q', action='store_true', help='"quiet" mode, disables automatic opening of graphs')
	parser.add_argument('-r', metavar='resolution', type=int, help='Resolution of graph in dpi. The default is 80dpi, which is suitable for viewing on monitors. Draft-quality images should be about 300dpi, and publication-quality images need to be 600 or 1200dpi depending on the journal.')
	parser.add_argument('-s', action='store_true', help='Force inputs to be interpreted as sequences (this is no longer a default behavior for infile args matching /[A-Z]+/')
	parser.add_argument('-t', metavar='format', default='png', help='Format of graph (\033[1mpng\033[0m, eps, jpeg, jpg, pdf, pgf, ps, raw, rgba, svg, svgz, tif, tiff')
	parser.add_argument('-v', action='store_true', help='Verbose output. Enables warnings and generally makes things messier')
	parser.add_argument('-w', '--walls', metavar='x(,dx(,y))', nargs='+', help='Draws bounds around sequences and such')
	parser.add_argument('-W', '--wall', metavar='x(,dx(,y))', nargs='+', help='Draws bounds around sequences and such')
	parser.add_argument('--axis-font', metavar='size', type=int, help='Axis label size (pt)')
	parser.add_argument('--height', metavar='height', type=float, help='Plot height in inches (default:5.5)')
	parser.add_argument('--mode', default='hydropathy', help='mode to run QUOD in (\033[1mhydropathy\033[0m, entropy)')
	parser.add_argument('--tick-font', metavar='size', type=int, help='Tick label size')
	parser.add_argument('--viewer', metavar='viewer', default=None, help='Viewer to be used for opening plots')
	parser.add_argument('--width', metavar='width', type=float, help='Plot width in inches (default:dynamic)')
	parser.add_argument('--window', metavar='windowsize', default=19, help='Window size for hydropathy')
	parser.add_argument('--x-offset', metavar='init_resi', default=0, type=int, help='Sets starting x-value')

	parser.add_argument('-ar', '--add-region', metavar='x0-x1(:color)(:"label")', nargs='+')

	parser.add_argument('-am', '--add-marker', metavar='(+id):x1,x2,x3,...xn(:color)', nargs='+', help='Adds a circular marker at the specified positions on the hydropathy curve of the specified sequence')

	parser.add_argument('-at', '--add-tms', metavar='x0-x1(:color)', nargs='+', help='Adds TMSs to plot, e.g. 40-60:red 65-90:blue to add a red TMS covering residues 40 to 60 and a blue one covering residues 65 to 90 (default color: orange)')
	parser.add_argument('-dt', '--delete-tms', metavar='(+id):x0-x1,x0-x1', nargs='+', help='Deletes TMSs on plot. Use +id to specify which sequence\'s TMSs should be deleted, e.g. "+0:5-25,30-40" to delete overlapping TMSs predicted for the first sequence and "+2:60-80:red,+4:70-95:blue" to delete overlapping TMSs predicted for the third and fifth sequence. +id specification propagates rightward and defaults to +0.')
	parser.add_argument('-et', '--extend-tms', metavar='(+id):x0-x1(:color)', nargs='+', help='Extends TMSs on plot. Use +id to specify which sequence\'s TMSs should be extended, e.g. "+0:5-25,30-40" to extend TMSs predicted for the first sequence to include residues 5-25 and 30-40 without changing colors and "+2:60-80:red,+4:70-95:blue" to extend TMSs predicted for the third sequence to include residues 60 to 80 and recolor them red and TMSs predicted for the fifth sequence to include residues 70 to 95 and recolor them blue. +id specification propagates rightward and defaults to +0.')
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

	main(args.infile, mode=args.mode, walls=args.walls, wall=args.wall, bars=args.bars, dpi=args.r, imgfmt=args.t, force_seq=args.s, outdir=args.d, outfile=args.o, color=args.color, title=args.title, quiet=args.q, viewer=args.a, axis_font=args.axis_font, width=args.width, height=args.height, x_offset=args.x_offset, add_tms=args.add_tms, delete_tms=args.delete_tms, extend_tms=args.extend_tms, replace_tms=args.replace_tms, no_tms=args.no_tms, tick_font=args.tick_font, add_marker=args.add_marker, add_region=args.add_region)

