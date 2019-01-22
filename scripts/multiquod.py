#!/usr/bin/env python

from __future__ import print_function, division, unicode_literals
import matplotlib
matplotlib.use('Qt4Agg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec
from collections import defaultdict
import shlex
import argparse
import os
import re
import numpy as np

import quod as quod

def blank_config():
	#TODO: encapsulate stuff to plot-specific classes
	cfg = {}
	cfg['n_rows'] = 0
	cfg['rowscheme'] = []
	cfg['subplots'] = set()
	cfg['sequences'] = defaultdict(lambda: '', {})
	cfg['mode'] = 'hvordan'
	cfg['outfile'] = None

	cfg['width'] = 8.5
	cfg['height'] = 11
	cfg['hres'] = 100
	cfg['hgap'] = 0.5
	cfg['vgap'] = -1.5
	#cfg['outdir'] = 'magnify_plots'
	cfg['dpi'] = 300
	cfg['weights'] = defaultdict(lambda: 1.0, {})
	cfg['color'] = {} #defaultdict(lambda: {'line':'red', 'tms':'orange'}))
	cfg['domains'] = {}
	cfg['label'] = {}
	cfg['linestyle'] = defaultdict(lambda: '-', {})
	cfg['xticks'] = {}
	cfg['title'] = {}

	cfg['walls'] = {}

	return cfg

def adj_font(args, cfg):
	name = args[0]
	if 'font' not in cfg: cfg['font'] = {}
	if name not in cfg['font']: cfg['font'][name] = {}
	try: cfg['font'][name][args[1]] = float(args[2])
	except ValueError: cfg['font'][name][args[1]] = args[2]

def set_prop(args, cfg):

	if len(args) < 1: raise TypeError('Not enough arguments')
	elif len(args) == 1:
		if False: pass
		else: raise TypeError('Unrecognized unary property {}'.format(args[0]))
	#figure-wide properties
	elif len(args) == 2:
		floatprops = (
			'width',
			'height',
			'hgap',
			'margin',
			'vgap',
		)
		intprops = (
			'hres',
			'dpi',
		)
		strprops = (
			'mode'
		)
		if args[0] in strprops: 
			cfg[args[0]] = args[1]
		elif args[0] in floatprops: 
			try: cfg[args[0]] = float(args[1])
			except ValueError as e: raise ValueError(e)
		elif args[0] in intprops: 
			try: cfg[args[0]] = int(args[1])
			except ValueError as e: raise ValueError(e)
		else: raise TypeError('Unrecognized property {}'.format(args[0]))
	#subplot-wide properties
	elif len(args) == 3:
		floatprops = (
			'xticks',
		)
		strprops = (
			'title',
			'xlabel',
			'ylabel',
		)

		if args[0] in strprops:
			if args[0] not in cfg: cfg[args[0]] = {}
			try: cfg[args[0]][args[1]] = args[2]
			except ValueError as e: raise ValueError(e)
		elif args[0] in floatprops:
			try: cfg[args[0]][args[1]] = float(args[2])
			except ValueError as e: raise ValueError(e)
	#entity-specific properties
	elif len(args) == 4:
		strprops = (
			'linestyle',
		)
		floatprops = (
			'linewidth',
		)
		if args[0] in strprops:
			try: cfg[args[0]][args[1]][int(args[2])] = args[3]
			except KeyError: cfg[args[0]][args[1]] = {int(args[2]):args[3]}
		elif args[0] in floatprops:
			if args[0] not in cfg: cfg[args[0]] = {}
			if args[1] not in cfg[args[0]]: cfg[args[0]][args[1]] = {}

			try: cfg[args[0]][args[1]][int(args[2])] = float(args[3])
			except KeyError: cfg[args[0]][args[1]] = {int(args[2]):float(args[3])}
		elif args[0] == 'font': adj_font(args[1:], cfg)
	elif len(args) >= 4:
		if args[0] == 'color':
			if args[1] not in cfg[args[0]]: cfg[args[0]][args[1]] = defaultdict(lambda: {'line':'red', 'tms':'orange'}, {})
			for colspec in args[3:]:
				if colspec.startswith('line:'): 
					cfg[args[0]][args[1]][int(args[2])]['line'] = colspec[5:]
				elif colspec.startswith('tms:'): 
					cfg[args[0]][args[1]][int(args[2])]['tms'] = colspec[4:]
				else: raise ValueError('Invalid color specification {}'.format(colspec))
		else: raise TypeError('Too many arguments')
	else: raise TypeError('Too many arguments')

	#elif len(args) < 2: raise TypeError('Unrecognized binary property {}'.format(args[0]))

	return cfg

def def_prop(args, cfg):
	if len(args) < 1: raise TypeError('Not enough arguments')
	else:
		name = args[0]
		cfg['subplots'].add(name)
		cfg['sequences'][name] = []
		cfg['label'][name] = name
		cfg['linestyle'][name] = defaultdict(lambda: '-', {})
		cfg['xticks']['name'] = None
		cfg['title'][name] = ''

		for i, arg in enumerate(args[1:]):
			if arg.startswith('seq:'):
				cfg['sequences'][name].append('>seq_{}\n{}'.format(i+1, arg[4:]))
			elif arg.startswith('file:'):
				with open(arg[5:]) as f:
					cfg['sequences'][name].append(f.read())
			else:
				with open(arg) as f:
					cfg['sequences'][name].append(f.read())
	return cfg

def color(args, cfg):
	name = args[0]
	seqi = int(args[1])
	for arg in args[2:]:
		if arg.startswith('line'): 
			cfg['color'][name][seqi]['line'] = arg[5:]
		elif arg.startswith('tms'): 
			cfg['color'][name][seqi]['tms'] = arg[4:]
	return cfg

def add_domain(args, cfg): 
	#, spans=[], yspan=[], label='', style='orange', alpha=None, pos='above', size=8, center=True)
	name = args[0]

	if len(args) < 6: raise IndexError('Incomplete domain specification (subplot, x1, x2, y1, y2, color[, label, valign, halign, size]')
	lims = [float(x) for x in args[1:5]]
	xlim = [lims[:2]]
	ylim = lims[2:4]
	ylim[1] -= ylim[0]
	color = args[5]

	#b 30 200 -2.8 -2.5 green "A domain" 8 above center
	#0  1   2    3    4     5         6  7     8      9
	if len(args) >= 7: label = args[6]
	else: label = 'unnamed_domain'

	if len(args) >= 8: font = float(args[7])
	else: font = 8.

	if len(args) >= 9: 
		valign = args[8]
		try: valign = float(valign)
		except ValueError: pass
	else: valign = 'above'

	if len(args) >= 10: halign = args[9]
	else: halign = 'center'
	#FIXME: add right-align support to quod
	if halign == 'center': halign = True
	else: halign = False

	kwargs = {'xlim':xlim, 'ylim':ylim, 'color':color, 'label':label, 'fontsize':font, 'pos':valign, 'halign':halign}

	try: cfg['domains'][name].append(kwargs)
	except KeyError: cfg['domains'][name] = [kwargs]
	return cfg
	
def add_walls(args, cfg):
	name = args[0]
	start = int(args[1])
	end = int(args[2])
	if len(args) >= 4: y = float(args[3])
	else: y = 2.

	if len(args) >= 5:
		if args[4] == '+': lim = 1
		elif args[4] == '-': lim = -1
		else: lim = None
	else: lim = None

	if len(args) >= 6: thickness = float(args[5])
	else: thickness = 1.0

	try: cfg['walls'][name].append({'start':start, 'end':end, 'y':y, 'lim':lim, 'thickness':thickness})
	except KeyError: cfg['walls'][name] = [{'start':start, 'end':end, 'y':y, 'lim':lim, 'thickness':thickness}]
	return cfg

def parse_config(f):
	cfg = blank_config()
	cfg['rawcfg'] = f.read()
	f.seek(0)

	for linenum, l in enumerate(f):
		if not l.strip(): continue
		elif l.strip().startswith('#'): continue
		cmd = shlex.split(l)
		if cmd[0] == 'addrow':
			cfg['rowscheme'].append([])
			cfg['n_rows'] += 1
			for name in cmd[1:]:
				cfg['rowscheme'][-1].append(name)
				cfg['subplots'].add(name)

		elif cmd[0] == 'set':
			try: set_prop(cmd[1:], cfg)
			except TypeError as e: raise TypeError('l. {}: set: {}'.format(linenum+1, e))
			except ValueError as e: raise ValueError('l. {}: set: {}'.format(linenum+1, e))

		elif cmd[0] == 'def':
			try: def_prop(cmd[1:], cfg)
			except TypeError as e: raise TypeError('l. {}: def: {}'.format(linenum+1, e))
			except ValueError as e: raise ValueError('l. {}: def: {}'.format(linenum+1, e))

		elif cmd[0] == 'addrow': cfg['rowscheme'].append(cmd[1:])

		#elif cmd[0] == 'color': 
		#	try: color(cmd[1:], cfg)
		#	except TypeError as e: raise TypeError('l. {}: color: {}'.format(linenum+1, e))
		#	except ValueError as e: raise ValueError('l. {}: color: {}'.format(linenum+1, e))

		elif cmd[0] == 'add':
			if cmd[1] == 'domain': add_domain(cmd[2:], cfg)
			elif cmd[1] == 'walls': add_walls(cmd[2:], cfg)
			else: raise TypeError('Unrecognized addition: {}'.format(cmd[1]))

		elif cmd[0] == 'tms':
			if 'tms' not in cfg: cfg['tms'] = []
			cfg['tms'].append(cmd[1:])

		elif cmd[0] == 'save':
			if len(cmd) < 2: raise TypeError('l. {}: save needs a filename'.format(linenum+1))
			cfg['outfile'] = cmd[1]

		else: raise TypeError('l. {}: Unrecognized directive: {}'.format(linenum+1, cmd[0]))

	#check that everything in rowscheme has been properly defined
	for row in cfg['rowscheme']:
		for name in row:
			if name not in cfg['subplots']:
				raise ValueError('Subplot {} not defined'.format(name))
	return cfg

def plot(cfg, name, fig, canvas, ax, entities, xlabel=' ', ylabel=None, yticks=None):
	quodplot = quod.Plot(fig=fig, canvas=canvas, ax=ax)
	quodplot.width = cfg['width']
	quodplot.height = cfg['height']
	for i, e in enumerate(entities):
		if type(e) is quod.What: 
			e.set_style(i)
			if name in cfg['color']:
				e.set_tms_color(cfg['color'][name][i]['tms'])
				e.set_curve_color(cfg['color'][name][i]['line'])
			if 'linestyle' in cfg and name in cfg['linestyle']:
				if i in cfg['linestyle'][name]: e.set_linestyle(cfg['linestyle'][name][i])
			if 'linewidth' in cfg and name in cfg['linewidth']:
				if i in cfg['linewidth'][name]: e.set_linewidth(cfg['linewidth'][name][i])
		quodplot.add(e)
	quodplot.render()

	ax.set_xlabel(xlabel)
	try: ax.set_title(cfg['titles'][name])
	except KeyError: ax.set_title('')
	if ylabel is not None: ax.set_ylabel(ylabel)
	if yticks is not None: ax.set_yticks(yticks)

	if name in cfg['xticks']:
		xlim = ax.get_xlim()
		ax.set_xticks(np.arange(xlim[0], xlim[1], cfg['xticks'][name]))

def do_tms_stuff(entities, cfg):

	if 'tms' not in cfg: return

	def list2line(l):
		s = ''
		for x in l: s += x + ' '
		return s.strip()
	for cmd in cfg['tms']:
		if cmd[0] == 'add':
			if len(cmd) < 5: raise TypeError('tms: Not enough arguments: {}'.format(list2line(cmd)))
			name = cmd[1]
			seqi = None
			color = None
			try: 
				seqi = int(cmd[2])
			except ValueError:
				color = cmd[2]
			start = 3 + (len(cmd[3:]) % 2)
			spans = []
			for i in range(start, len(cmd), 2):
				spans.append([float(x) for x in cmd[i:i+2]])
			if seqi is None:
				entities[name].append(quod.HMMTOP('', nohmmtop=True, style=color))
				entities[name][-1].spans = spans
			else:
				si = 0
				found = False
				for ei, e in enumerate(entities[name]):
					if si == seqi:
						e.entities[1].spans += spans
						found = True
						break
					if type(e) is quod.What: si += 1
				if not found: 
					raise ValueError('tms: Could not find subplot {} sequence {}: {}'.format(name, seqi, list2line(cmd)))

		elif cmd[0] in ('load', 'append'):
			if len(cmd) < 4: raise TypeError('tms: Not enough arguments: {}'.format(list2line(cmd)))
			name = cmd[1]
			seqi, color = None, None
			try: seqi = int(cmd[2])
			except ValueError: color = cmd[2]
			spans = []
			for rawfn in cmd[3:]:
				if rawfn.startswith('file:'): fn = rawfn[5:]
				else: fn = rawfn
				with open(fn) as f:
					for l in f:
						raw = re.findall('(?:[0-9]+\s*)+$', l.strip())
						if not raw: continue
						else:
							for s in raw: 
								indices = [int(x) for x in s.split()]
								for i in range(len(indices)%2, len(indices), 2):
									spans.append(indices[i:i+2])
						#print(re.findall('(?:[0-9]+\s*)+$', l.strip()))
			if seqi is None:
				entities[name].append(HMMTOP('', style=color, nohmmtop=True))
				entities[name][-1].spans = spans
			else:
				si = 0
				found = False
				for ei, e in enumerate(entities[name]):
					if si == seqi:
						if cmd[0] == 'load': e.entities[1].spans = spans
						elif cmd[0] == 'append': e.entities[1].spans += spans
						else: raise Exception('Impossible exception')
						found = True
						break
					if type(e) is quod.What: si += 1
				if not found: 
					raise ValueError('tms: Could not find subplot {} sequence {}: {}'.format(name, seqi, list2line(cmd)))

		elif cmd[0] == 'erase':
			if len(cmd) < 5: raise TypeError('tms: Not enough arguments: {}'.format(list2line(cmd)))
			name = cmd[1]
			spans = []

			for i in range(3, len(cmd), 2): spans.append([int(x) for x in cmd[i:i+2]])

			si = 0
			seqi = int(cmd[2])
			found = False

			#TODO: automerge eraser targets

			for ei, e in enumerate(entities[name]):
				if si == seqi:

					popme = []
					appendme = []
					for i, oldspan in enumerate(e.entities[1].spans):
						for eraseme in spans:
							#case 1: oldspan is a subset of eraseme
							if (eraseme[0] <= oldspan[0] and oldspan[1] <= eraseme[1]):
								popme.append(i)
							#case 2: eraseme is a subset of oldspan
							elif (eraseme[0] >= oldspan[0] and oldspan[1] >= eraseme[1]):
								popme.append(i)
								appendme.append([oldspan[0], eraseme[0]])
								appendme.append([oldspan[1], eraseme[1]])
							#case 3: eraseme overlaps oldspan
							elif (eraseme[0] <= oldspan[0] <= eraseme[1]): 
								e.entities[1].spans[i][0] = max(oldspan[0], eraseme[1])
							#case 4: oldspan overlaps eraseme 
							if (eraseme[0] <= oldspan[1] <= eraseme[1]): 
								e.entities[1].spans[i][0] = min(oldspan[1], eraseme[0])
					popme = sorted(set(popme))
					for i in popme[::-1]: e.entities[1].spans.pop(i)

					e.entities[1].spans += appendme
					found = True
				if type(e) is quod.What: si += 1
			if not found: 
				raise ValueError('tms: Could not find subplot {} sequence {}: {}'.format(name, seqi, list2line(cmd)))

		elif cmd[0] == 'delete':
			if len(cmd) < 4: raise TypeError('tms: Not enough arguments: {}'.format(list2line(cmd)))
			name = cmd[1]
			seqi = int(cmd[2])
			si = 0
			found = False

			deleteme = sorted([int(x) for x in cmd[3:]])[::-1]
			for ei, e in enumerate(entities[name]):
				if si == seqi:
					for i in deleteme: e.entities[1].spans.pop(i)
					found = True
					break
				if type(e) is quod.What: si += 1
			if not found: 
				raise ValueError('tms: Could not find subplot {} sequence {}: {}'.format(name, seqi, list2line(cmd)))


	pass

def run_quod(cfg):

	entities = {}
	maxlens = {}
	for name in cfg['subplots']:
		entities[name] = []
		maxlens[name] = 0
		for seq in cfg['sequences'][name]:
			entities[name].append(quod.What(seq))
			maxlens[name] = max(maxlens[name], len(entities[name][-1]))


			#walls
			if name in cfg['walls']:
				for walldef in cfg['walls'][name]:
					entities[name].append(quod.Wall(spans=[[walldef['start'], walldef['end']]], y=walldef['y'], ylim=walldef['lim'], thickness=walldef['thickness'], wedge=walldef['thickness']))

			#domains
			if name in cfg['domains']:
				for domaindef in cfg['domains'][name]:
					entities[name].append(quod.Region(
						spans=domaindef['xlim'], 
						yspan=domaindef['ylim'],
						label=domaindef['label'],
						style=domaindef['color'],
						pos=domaindef['pos'],
						size=domaindef['fontsize'],
						center=domaindef['halign']
					))

		cfg['weights'][name] *= maxlens[name]

	
	if 'tms' in cfg: do_tms_stuff(entities, cfg)		

	fig = Figure()
	canvas = FigureCanvas(fig)
	fig.set_tight_layout(True)

	gs = gridspec.GridSpec(len(cfg['rowscheme']), cfg['hres'])
	hgap = cfg['hgap']/cfg['width']/2
	margin = 0.03 if 'margin' not in cfg else cfg['margin']
	gs.update(
		left=margin + hgap, 
		right=cfg['width'] - margin, 
		top=cfg['height'] - margin-hgap, 
		bottom=margin, wspace=0)

	def get_limit(wt1, wtlist, offset=0): 
		x = int((offset + wt1/sum(wtlist)) * cfg['hres'])
		return max(0, min(x, cfg['hres']))

	if cfg['mode'] == 'hvordan':

		lims = []
		cfg['lims'] = {}
		for namelist in cfg['rowscheme']:
			row = [cfg['weights'][name] for name in namelist]
			last = -hgap
			if len(namelist) == 1: 
				cfg['lims'][namelist[0]] = [0, cfg['hres']]
			else:
				bounds = [0]

				s = 0
				rawweights = [cfg['weights'][name] for name in namelist]
				weights = [w/sum(rawweights) for w in rawweights]

				for i in range(len(namelist) - 1):
					bounds.append((weights[i] - hgap) * cfg['hres'])
					bounds.append((weights[i] + hgap) * cfg['hres'])

				#for i in namelist


				bounds.append(cfg['hres'])
				bounds = [int(x) for x in bounds]
				
				for i in range(0, len(namelist)):
					name = namelist[i]
					cfg['lims'][name] = [bounds[2*i], bounds[2*i + 1]]



		axdict = {}
		for r, row in enumerate(cfg['rowscheme']):
			for name in row:
				#axdict[name] = fig.add_subplot(gs[r, cfg['lims'][name][0]:cfg['lims'][name][1]])
				axdict[name] = fig.add_subplot(gs[r,cfg['lims'][name][0]:cfg['lims'][name][1]])
			

		n = 0
		name2row = {}
		stretchlabel = {}
		firstcol = set()

		for r, row in enumerate(cfg['rowscheme']):
			for c, name in enumerate(row):
				if not c: firstcol.add(name)
				name2row[name] = r
				if 'ylabel' in cfg and name in cfg['ylabel']: ylabel = cfg['ylabel'][name]
				else: ylabel = '' if c else None
				yticks = [] if c else None
				plot(cfg, name, fig, canvas, axdict[name], entities[name], ylabel=ylabel, yticks=yticks)
				try: title = cfg['title'][name]
				except KeyError: title = chr(65 + n)
				axdict[name].set_title(title, loc='left')
				n += 1

			if len(row) > 1:
				ax = fig.add_subplot(gs[r,:])
				ax.set_xticks([0])
				ax.axes.get_yaxis().set_visible(False)
				ax.set_frame_on(False)
				#TODO: expose customizing xlabels per-plot/per-row
				ax.set_xlabel('Position')
				ax.tick_params(labelcolor=(0,0,0,0), color=(0,0,0,0))
				stretchlabel[r] = ax
			elif len(row) == 1:
				#ax = fig.add_subplot(gs[r,:])
				ax = axdict[row[0]]
				ax.set_xlabel('Position')
				stretchlabel[r] = ax

		if 'xlabel' in cfg:
			for name in cfg['xlabel']:
				stretchlabel[name2row[name]].set_xlabel(cfg['xlabel'][name])
		if 'ylabel' in cfg:
			for name in cfg['ylabel']:
				stretchlabel[name2row[name]].set_ylabel(cfg['ylabel'][name])

		if 'font' in cfg:
			for name in cfg['font']:
				for target in cfg['font'][name]:
					if target.endswith('ticks'):
						axdict[name].tick_params(labelsize=cfg['font'][name][target])
						if name in firstcol:
							stretchlabel[name2row[name]].tick_params(labelsize=cfg['font'][name][target])
					elif target.endswith('xaxis'):
						axdict[name].set_xlabel(axdict[name].get_xlabel(), 
							fontsize=cfg['font'][name][target])
						stretchlabel[name2row[name]].set_xlabel(stretchlabel[name2row[name]].get_xlabel(), 
							fontsize=cfg['font'][name][target])
					elif target.endswith('yaxis'):
						axdict[name].set_ylabel(axdict[name].get_ylabel(), 
							fontsize=cfg['font'][name][target])
					elif target == 'title':
						axdict[name].set_title(axdict[name].get_title(),
							fontsize=cfg['font'][name][target])


		#gs.tight_layout(fig, pad=cfg['margin'], w_pad=cfg['hgap'], h_pad=0.0)
		gs.tight_layout(fig, pad=0, w_pad=cfg['hgap'], h_pad=cfg['vgap'])
		fig.savefig(cfg['outfile'], dpi=cfg['dpi'])

	else: raise NotImplementedError('Unimplemented mode: {}'.format(cfg['mode']))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', help='configuration file to load')
	args = parser.parse_args()

	with open(args.infile) as f:
		cfg = parse_config(f)
	if not cfg['outfile']: raise ValueError('No outfile defined ("save FILENAME")')
	run_quod(cfg)


