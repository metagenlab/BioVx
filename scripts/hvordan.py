#!/usr/bin/env python2

from __future__ import print_function
import os, sys, subprocess, re 

import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import time, hashlib
import tempfile

import numpy as np

import quod, tcblast

#import Bio.Entrez
import Bio.pairwise2, Bio.SubsMat.MatrixInfo

DEBUG = 0
VERBOSITY = 1

def warn(*msgs):
	for l in msgs: print('[WARNING]', l, file=sys.stderr)
def error(*msgs):
	for l in msgs: print('[ERROR]', l, file=sys.stderr)
	exit()
def info(*msgs):
	for l in msgs: print('[INFO]', l, file=sys.stderr)

def fetch(accessions, email=None, db='protein'):
	if not accessions: return ''
	if db == 'tcdb':
		out = ''
		for acc in accessions:
			try: 
				if DEBUG: info('Running blastdbcmd')
				fa = subprocess.check_output(['blastdbcmd', '-db', 'tcdb', '-target_only', '-entry', acc])
				out += fa + '\n'
			except ValueError: raise ValueError
		return out
	else:
		if DEBUG: info('Preparing to fetch non-TCDB sequences')
		acclist = ''
		for x in accessions: acclist += ',' + x
		acclist = acclist[1:]

		try:
			if DEBUG: info('Running blastdbcmd')
			p = subprocess.Popen(['blastdbcmd', '-db', 'nr', '-entry', acclist], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			#out = re.sub('>', '\n', out) + '\n'

			remotes = ''
			for l in err.split('\n'):
				if l.strip(): remotes += '%s,' % l.split()[-1]
			remotes = remotes[:-1]
			if DEBUG: info('Fetching from remote')
			#out += subprocess.check_output(['curl', '-d', 'db=%s&id=%s&rettype=fasta&retmode=text' % (db, acclist), 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'])
			if remotes: out += subprocess.check_output(['curl', '-d', 'db=%s&id=%s&rettype=fasta&retmode=text' % (db, remotes), 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'])
			#out += subprocess.check_output(['curl', '-d', 'db=protein&id=Q9RBJ2&rettype=fasta&retmode=text', 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'])

			return out
		except subprocess.CalledProcessError:

			if not email:
				if 'ENTREZ_EMAIL' in os.environ: email = os.environ['ENTREZ_EMAIL']
				else: 
					raise TypeError('Missing argument email')

			if DEBUG: info('Fetching from remote')
			out += subprocess.check_output(['curl', '-d', 'db=%s&id=%s&rettype=fasta&retmode=text' % (db, acclist), 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'])
			return out
	if DEBUG: info('Done fetching a batch from %s' % db)

def parse_p2report(p2report, minz=15, maxz=None, musthave=None, thispair=None):

	if musthave and thispair: error('Arguments musthave and thispair are not mutually compatible')

	line = 0

	if minz == None: minz = -2**16
	if maxz == None: maxz = 2**16

	bcs = []
	alnregs = {}
	stats = {}
	for l in p2report.split('\n'):
		line += 1
		if line == 1:
			fams = [re.split('[, :]+', l)[2], re.split('[, :]+', l)[4]]
			bcs = {fams[0]:[], fams[1]:[]}
		elif line == 2: pass
		else: 
			if not l.strip(): continue
			ls = l.split('\t')
			z = float(ls[3])
			if minz <= z <= maxz: 
				#bcs.append(ls[:2])
				if musthave and ls[0] not in musthave and ls[1] not in musthave: continue
				found = 1
				if thispair:
					found = 0
					for pair in thispair:
						if ls[:2] != pair and ls[:2][::-1] != pair: continue
						else: 
							found = 1
							break
				if not found: continue
				bcs[fams[0]].append(ls[0])
				bcs[fams[1]].append(ls[1])
				try: alnregs[ls[0]][ls[1]] = (ls[6], ls[7])
				except KeyError: alnregs[ls[0]] = {ls[1]:(ls[6], ls[7])}
				try: stats[ls[0]][ls[1]] = ls[2:6]
				except KeyError: stats[ls[0]] = {ls[1]:ls[2:6]}
			
	return fams, bcs, alnregs, stats

def seek_initial(p1ds, bcs):
	hits = {}
	for fam in sorted(bcs):
		hits[fam] = {}
		for bc in sorted(bcs[fam]): hits[fam][bc] = []

	fs = {}
	for p1d in p1ds:
		if os.path.isfile(p1d):
			for bc in sorted(bcs): fs[bc] = p1d
		elif os.path.isdir(p1d):
			for fam in sorted(bcs):
				if os.path.isfile('%s/%s.tbl' % (p1d, fam)): fs[fam] = '%s/%s.tbl' % (p1d, fam)

				elif os.path.isfile('%s/%s/psiblast.tbl' % (p1d, fam)): fs[fam] = '%s/%s/psiblast.tbl' % (p1d, fam)

				else: error('Could not find famXpander results table in %s' % p1d)

		else: error('Could not find p1d %s' % p1d)
	for bc in sorted(bcs):
		with open(fs[bc]) as f:
			for l in f:
				if not l.strip(): continue
				if l.lstrip().startswith('#'): continue
				if '\t' not in l: continue
				ls = l.split('\t')

				try: hits[bc][ls[1]].append((float(ls[4]), ls[0], (int(ls[6]), int(ls[7])), (int(ls[9]), int(ls[10]))))
				except KeyError: hits[bc][ls[1]] = [(float(ls[4]), ls[0], (int(ls[6]), int(ls[7])), (int(ls[9]), int(ls[10])))]

	for fam in sorted(bcs):
		for bc in sorted(hits[fam]): hits[fam][bc] = sorted(hits[fam][bc])[0]
	return hits

def clean_fetch(accs, outdir, force=False, email=None):
	if DEBUG: info('Fetching %s' % accs)
	if not force:
		removeme = []
		for acc in accs:
			if os.path.isfile(outdir + '/%s.fa' % acc): removeme.append(acc)
		for acc in removeme: accs.remove(acc)

	if not os.path.isdir(outdir): os.mkdir(outdir)

	dlme = []
	tcdlme = []
	for acc in accs:
		if os.path.isfile(outdir + '/%s.fa' % acc): continue
		else: 
			if re.match('[0-9]\.[A-Z]\.[0-9]+\.', acc): tcdlme.append(acc)
			else: dlme.append(acc)

	allfaa = ''
	if dlme: 
		if VERBOSITY: info('Downloading %d sequence(s)' % len(dlme))
		allfaa += fetch(dlme, email=email)

	if tcdlme:
		if VERBOSITY: info('Loading %d TCDB sequence(s)' % len(tcdlme))
		allfaa += fetch(tcdlme, db='tcdb', email=email)

		if VERBOSITY: info('Done loading %d TCDB sequence(s)' % len(tcdlme))

	with open('%s/allseqs.faa' % outdir, 'w') as f: f.write(allfaa)

	with open('%s/allseqs.faa' % outdir) as f: 
		faa = Bio.SeqIO.parse(f, format='fasta')
		for record in faa: 
			for desc in record.description.split('>'):
				name = desc.split()[0]
				if name.count('.') < 4: name = name[:name.find('.')]
				if name.count('|') == 1: name = name.split('|')[1]

				if DEBUG > 1: info('Saving %s' % name)
				with open('%s/%s.fa' % (outdir, name), 'w') as f: f.write('>%s\n%s' % (desc, record.seq))

	fastas = {}

	#for fa in allfaa.split('\n\n'):
	#	if not fa.strip(): continue
	#	for acc in accs:
	#		if acc in fastas: pass
	#		if fa.startswith('>' + acc):
	#			fastas[acc] = fa

	#for x in sorted(fastas): 
	#	if DEBUG: info('Saving %s' % x)
	#	f = open(outdir + '/%s.fa' % x, 'w')
	#	f.write(fastas[x])
	#	f.close()

def quod_set(seqids, sequences, indir, outdir, dpi=300, force=False, bars=[], prefix='', suffix='', silent=False, pars=[]):
	if not os.path.isdir(outdir): os.mkdir(outdir)

	#wedges = [[[x, 2 * (0.5 - (i % 2))] for i, x in enumerate(span)] for span in bars]

	ove = lambda x: int(2 * (0.5 - (x % 2)))

	wedges = []
	for i, span in enumerate(bars):
		wedges.append([])
		for j, x in enumerate(span):
			#wedges.append(quod.Wedge(
			if 1 <= i <= 2: y = -2
			else: y = 0
			wedges[-1].append(quod.Wedge(x=x, dx=ove(j), y=y))

	medges = []
	for i, span in enumerate(pars):
		medges.append([])
		for j, x in enumerate(span):
			y = 2
			medges[-1].append(quod.Wedge(x=x, dx=ove(j), y=y))

	#Draw A: barred by B
	quod.what([sequences[seqids[0]]], title=seqids[0], imgfmt='png', directory=outdir, filename=(seqids[0] + '_' + seqids[1] + '.png'), dpi=dpi, hide=1, bars=bars[0], wedges=wedges[0], silent=silent)

	#Draw B: barred by C
	quod.what([sequences[seqids[1]]], title=seqids[1], imgfmt='png', directory=outdir, filename=(seqids[1] + '_' + seqids[2] + '.png'), dpi=dpi, hide=1, bars=bars[1], wedges=wedges[1]+medges[0], silent=True)

	#Draw C: barred by B
	quod.what([sequences[seqids[2]]], title=seqids[2], imgfmt='png', directory=outdir, filename=(seqids[2] + '_' + seqids[1] + '.png'), dpi=dpi, hide=1, bars=bars[2], color=1, wedges=wedges[2]+medges[1], silent=True)

	#Draw D: barred by C
	quod.what([sequences[seqids[3]]], title=seqids[3], imgfmt='png', directory=outdir, filename=(seqids[3] + '_' + seqids[2] + '.png'), dpi=dpi, hide=1, bars=bars[3], color=1, wedges=wedges[3], silent=True)


def build_html(bc, indir, blasts, outdir='hvordan_out/html', filename='test.html', lastpair=None, nextpair=None):

	if not os.path.isdir(outdir): os.mkdir(outdir)
	if not os.path.isdir(outdir + '/assets'): os.mkdir(outdir + '/assets')

	if not os.path.isfile(outdir + '/assets/openclose.js'):
		f = open(outdir + '/assets/openclose.js', 'w')
		f.write('function toggle_section(sectionid, selfid) {\n\tvar section = document.getElementById(sectionid);\n\tvar me = document.getElementById(selfid);\n\t//console.log([section, section.style.display]);\n\tif (section.style.display == \'none\') {\n\t\tsection.style.display = \'block\';\n\t\tme.innerHTML = \'Hide\';\n\t} else { \n\t\tsection.style.display = \'none\'; \n\t\tme.innerHTML = \'Show\';\n\t}\n}')
		f.close()
	if not os.path.isfile(outdir + '/assets/nice.css'):
		f = open(outdir + '/assets/nice.css', 'w')
		f.write('body {\n\tfont-family: sans-serif;\n\theight: 100%;\n}\ndiv {\n\tdisplay: block;\n}\ndiv.tcblast {\n\tmax-width: 1500px;\n}\ndiv.fullblast {\n\twidth: 50%;\n\tfloat: left;\n}\ndiv.tabular1 {\n\twidth: 49%;\n\tfloat: left;\n\theight: 100%;\n}\ndiv.tabular2 {\n\twidth: 49%;\n\tfloat: right;\n\theight: 100%;\n}\nimg.bluebarplot {\n\tmax-width: 100%;\n\theight: auto;\n}\n.clear { clear: both; }\n.scrollable {\n\toverflow-y: scroll;\n}\n.resizeable {\n\tresize: vertical;\n\toverflow: auto;\n\tborder: 1px solid gray;\n\tdisplay: block;\n\tpadding-bottom: 1ex;\n}\n.bluebars {\n\theight: 25vh;\n}\n.pairwise {\n\theight: 50vh;\n}\n.whatall {\n\theight: 50vh;\n}\n.whataln {\n\twidth: 100%;\n}\n#seqs {\n\tdisplay: none;\n}\n\n\n\n.summtbl {\n\tfont-family: monospace, courier;\n\tfont-size: 75%;\n}\n.oddrow {\n\tbackground-color: #d8d8d8;\n}\ntd {\n\tpadding-right: 1em;\n}\n.red {\n\tcolor: red;\n}\nimg {\n\tborder: 1pt solid black;\n}\n.monospace {\n\tfont-family: monospace;\n}')
		f.close()
	#bc := [WP_1234567890, AP_1234567890]
	title = 'HVORDAN summary: %s vs %s' % tuple(bc[1:3])

	out = '<html><head><title>%s</title>' % title
	out += '\n<link rel="stylesheet" type="text/css" href="assets/nice.css"/>'
	out += '\n<script src="assets/openclose.js"></script>'
	out += '\n</head><body>'

	out += '\n<h1>%s</h1>' % title
	if lastpair or nextpair:
		out += '\n'
		if lastpair: out += '<a href="%s_vs_%s.html">&#9664; %s vs %s</a> ' % (lastpair[1], lastpair[2], lastpair[1], lastpair[2])
		if nextpair: out += '<a href="%s_vs_%s.html">%s vs %s &#9654;</a> ' % (nextpair[1], nextpair[2], nextpair[1], nextpair[2])
		out += '<br/>'
	out += '\n<h2>Table of contents</h2>'
	out += '\n<button class="showhide" id="tocsh" onclick="toggle_section(\'toc\', \'tocsh\')">Hide</button>'

	out += '\n<div class="toc" id="toc"> <ol> <li><a href="#summary">Summary</a></li> <li><a href="#tcsummary">TCBLAST Summary</a></li> <li><a href="#pairwise">Pairwise</a></li> <li><a href="#abcd">ABCD hydropathy plots</a></li> <li><a href="#bc">BC hydropathy plot</a></li> <li><a href="sequences">Sequences</a></li> </ol> </div>'

	#stats
	out += '\n<h2>Summary</h2>'
	out += '\n<button class="showhide" id="summarysh" onclick="toggle_section(\'summary\', \'summarysh\')">Hide</button>'
	out += '\n<div class="whataln" id="summary">'
	out += '\nSS Z-score: %s<br/>' % bc[8]
	out += '\nGSAT Z-score: %s<br/>' % bc[9]
	out += '\nSubject align-length: %s<br/>' % bc[10]
	out += '\nTarget align-length: %s<br/>' % bc[11]
	out += '\n</div>'

	out += '\n<h2>TCBLAST</h2>'

	#bluebars
	out += '\n<button class="showhide" id="tcblastsh" onclick="toggle_section(\'tcblast\', \'tcblastsh\')">Hide</button>'
	out += '\n<div class="tcblast" id="tcblast"><a name="tcsummary"><h3>TCBLAST Summary</h3></a>'
	out += '\n<div class="resizeable bluebars"><div class="scrollable tabular1">'
	out += '\n<img class="bluebarplot" src="../graphs/TCBLAST_%s.png"/>' % bc[1]
	out += '\n</div><div class="scrollable tabular2">'
	out += '\n<img class="bluebarplot" src="../graphs/TCBLAST_%s.png"/>' % bc[2]
	out += '\n</div></div>'

	#pairwise
	out += '\n<div class="clear"></div><a name="pairwise"><h3>Pairwise</h3></a><div class="resizeable pairwise"><div class="scrollable tabular1">'
	out += '\n%s' % blasts[0][1]
	out += '</div><div class="scrollable tabular2">'
	out += '\n%s' % blasts[1][1]
	out += '\n</div></div></div>'


	#abcd bc
	out += '\n<div class="clear"></div><a name="abcd"><h3>ABCD Hydropathy plots</h3></a>'
	out += '\n<button class="showhide" id="abcdsh" onclick="toggle_section(\'abcd\', \'abcdsh\')">Hide</button>'

	out += '\n<div class="whatall" id="abcd">'
	out += '\n<div class="tabular1">'
	out += '\nA<br/><img class="bluebarplot" id="plota" src="../graphs/%s_%s.png"/><br/>' % (bc[0], bc[1])
	out += '\nB<br/><img class="bluebarplot" id="plotb" src="../graphs/%s_%s.png"/><br/>' % (bc[1], bc[2])
	out += '\n</div><div class="tabular2">'
	out += '\nD<br/><img class="bluebarplot" id="plotd" src="../graphs/%s_%s.png"/><br/>' % (bc[3], bc[2])
	out += '\nC<br/><img class="bluebarplot" id="plotc" src="../graphs/%s_%s.png"/><br/>' % (bc[2], bc[1])
	out += '\n</div></div>'

	out += '\n<div class="clear"></div><br/><a name="bc"><h3>BC hydropathy plot</h3></a>'
	out += '\n<button class="showhide" id="bcsh" onclick="toggle_section(\'bc\', \'bcsh\')">Hide</button>'
	out += '\n<div class="resizeable whataln" id="bc"><div class="scrollable">'
	out += '<img class="bluebarplot" id="plotbc" src="../graphs/%s_vs_%s.png"/><br/>' % (bc[1], bc[2])
	out += '\n</div></div>'

	#out += '\n<button class="showhide" id="tcblastsh" onclick="toggle_section(\'tcblast\', \'tcblastsh\')">Hide</button>'

	out += '\n<br/><div style="height: 10ex"></div>'

	#sequences
	out += '\n<div class="clear"></div><br/><a name="sequences"><h3>Sequences</h3></a>'
	out += '\n<button class="showhide" id="seqsh" onclick="toggle_section(\'sequences\', \'seqsh\')">Hide</button>'
	out += '\n<div class="resizeable whataln monospace" id="sequences"><div class="scrollable">'
	out += ('\n%s\n%s\n%s\n%s' % tuple(bc[4:8])).replace('\n', '<br/>\n')
	out += '\n</div></div>'

	out += '\n</body></html>'

	f = open(outdir + '/' + filename, 'w')
	f.write(out)
	f.close()

def get_fulltrans(fams, bcs, abcd):

	pairs = zip(bcs[fams[0]], bcs[fams[1]])
	origs = [abcd[fams[0]], abcd[fams[1]]]

	fulltrans = []
	for p in pairs:
		fulltrans.append(tuple([origs[0][p[0]][1], p[0], p[1], origs[1][p[1]][1]]))
	return fulltrans

def blastem(acc, indir, outdir, dpi=300, force=False, seqbank={}, tmcount={}, maxhits=50):
	f = open(indir + '/sequences/' + acc + '.fa')
	seq= f.read()
	f.close()

	return tcblast.til_warum(seq, outfile='%s/graphs/TCBLAST_%s.png' % (outdir, acc), title=acc, dpi=dpi, outdir='%s/blasts' % outdir, clobber=force, seqbank=seqbank, tmcount=tmcount, silent=True, maxhits=maxhits)

	#fn = outdir + '/' + filename + '.png'

	#blasts = tcblast.til_warum(seq, fn, dpi=dpi)
	#blasts = [tcblast.til_warum(l[0], args.o + '/images/' + accs[0] + '.png', dpi=args.r, html=2, outdir=args.o + '/hmmtop'), tcblast.til_warum(l[1], args.o + '/images/' + accs[1] + '.png', dpi=args.r, html=2, outdir=args.o + '/hmmtop')]

def identifind(seq1, seq2):
	#Seq1 = Bio.Seq.Seq(seq1, Bio.Alphabet.ProteinAlphabet())
	if seq1.startswith('>'): seq1 = seq1[seq1.find('\n')+1:]
	if seq2.startswith('>'): seq2 = seq2[seq2.find('\n')+1:]
	seq1 = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', seq1.upper())
	seq2 = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', seq2.upper())

	if DEBUG: info('Starting an alignment')
	#alns = Bio.pairwise2.align.localds(seq1, seq2, Bio.SubsMat.MatrixInfo.ident, -10, -0.5)
	#out = subprocess.check_output(['ggsearch36'])
	aln = ggsearch(seq1, seq2)

	if DEBUG: info('Finished an alignment')

	subjstart = 0

	#sngap = re.findall('^-+', aln[0])
	#if sngap: sngap = len(sngap[0])
	#else: sngap = 0

	#scgap = re.findall('-+$', aln[0])
	#if scgap: scgap = len(aln[0]) - len(scgap[0]) - 1
	#else: scgap = len(aln[0])-1

	#tngap = re.findall('^-+', aln[1])
	#if tngap: tngap = len(tngap[0])
	#else: tngap = 0

	#tcgap = re.findall('-+$', aln[1])
	#if tcgap: tcgap = len(aln[1]) - len(tcgap[0]) - 1
	#else: tcgap = len(aln[1])-1

	#if sngap: 
	#	sstart = 0
	#	tstart = sngap
	#else: 
	#	sstart = tngap
	#	tstart = 0

	igap1 = re.findall('^-+', aln[0])
	igap2 = re.findall('^-+', aln[1])
	tgap1 = re.findall('-+$', aln[0])
	tgap2 = re.findall('-+$', aln[1])
	if igap1:
		#1 -----CYFQNCPRG
		#2 CYFQNCPRGCYFQN
		qstart = 0
		sstart = len(igap1[0])
	elif igap2:
		#1 CYFQNCPRGCYFQN
		#2 -----CYFQNCPRG
		qstart = len(igap2[0])
		sstart = 0
	else:
		#1 CYFQNCPRGCYFQN
		#2 CYFQNCPRG-----
		qstart = 0
		sstart = 0
	if tgap1:
		#1 CYFQNCPRG-----
		#2 CYFQNCPRGCYFQN
		qend = len(seq1)-1
		send = len(seq2)-1-len(tgap1[0])
	elif tgap2:
		#1 CYFQNCPRGCYFQN
		#2 CYFQNCPRG-----
		qend = len(seq1)-1-len(tgap[1])
		send = len(seq2)-1
	else:
		#1 CYFQNCPRGCYFQN
		#2 -----CYFQNCPRG
		qend = len(seq1)-1
		send = len(seq2)-1

	return qstart+1, qend+1, sstart+1, send+1

		#I prefer 0-indexing, but pretty much everyone 1-indexes (at least for protein sequences)

def ggsearch(seq1, seq2):
	if not seq1.startswith('>'): seq1 = '>seq1\n' + seq1
	if not seq2.startswith('>'): seq2 = '>seq2\n' + seq2

	try:
		f1 = tempfile.NamedTemporaryFile(delete=False)
		f1.write(seq1)
		f1.close()

		f2 = tempfile.NamedTemporaryFile(delete=False)
		f2.write(seq2)
		f2.close()

		cmd = ['ggsearch36', '-m', '3', f1.name, f2.name]
		out = subprocess.check_output(cmd)

	finally:
		os.remove(f1.name)
		os.remove(f2.name)
	
	seqi = 0
	alns = []
	for l in out.split('\n'):
		if l.startswith('>'): seqi += 1
		if seqi:
			if not l.strip(): seqi = 0
			#elif l.startswith('>'): alns.append(l + '\n')
			#else: alns[-1] += l + '\n'
			elif l.startswith('>'): alns.append('')
			else: alns[-1] += l

	return alns


def summarize(p1d, p2d, outdir, minz=15, maxz=None, dpi=100, force=False, email=None, musthave=None, thispair=None, fams=None, maxhits=50):
	if thispair is not None:
		if len(thispair) % 2: error('Unpaired sequence found')
		else:
			truepairs = [thispair[i:i+2] for i in range(0, len(thispair), 2)]
	else: truepairs = None

	if not os.path.isdir(outdir): os.mkdir(outdir)

	if VERBOSITY: info('Reading Protocol2 report')
	try: f = open(p2d + '/report.tbl')
	except IOError:
		if os.path.isfile(p2d):
			f = open(p2d)
			warn('Opening %s as a Protocol2 results table' % p2d)
		else:
			try:
				famvfam = '%s_vs_%s' % tuple(fams)
				try: 
					f = open('%s/%s/report.tbl' % (p2d, famvfam))
					info('Could not find report.tbl in %s, falling back on family vs family subdirectory' % p2d)
				except IOError:
					try: f = open('%s/%s/%s/report.tbl' % (p2d, famvfam, famvfam))
					except IOError: error('Could not find a Protocol2 directory for %s and %s' % tuple(fams))
			except TypeError: error('Specify families if using Protocol2 root directories')
	p2report = f.read()
	f.close()

	fams, bcs, alnregs, stats = parse_p2report(p2report, minz, maxz, musthave=musthave, thispair=truepairs)

	if VERBOSITY: info('Selecting best A-B C-D pairs')
	abcd = seek_initial(p1d, bcs)
	#for k in abcd: 
	#	for j in abcd[k]: 
	#		print(k, j, abcd[k][j])
	fulltrans = get_fulltrans(fams, bcs, abcd)

	fetchme = set()
	pairstats = {}
	#for fam in abcd:
	#	for bc in abcd[fam]:
	#		fetchme.add(bc) # B|C
	#		fetchme.add(abcd[fam][bc][1]) #A|D
	#		try: pairstats[bc][abcd[fam][bc][1]] = abcd[fam][bc]
	#		except KeyError: pairstats[bc] = {abcd[fam][bc][1]:abcd[fam][bc]}
	for fam in abcd:
		for bc in abcd[fam]:
			try: pairstats[bc][abcd[fam][bc][1]] = abcd[fam][bc]
			except KeyError: pairstats[bc] = {abcd[fam][bc][1]:abcd[fam][bc]}

	for pair in fulltrans:
		for acc in pair: fetchme.add(acc)

	#grab all relevant sequences and store them
	if VERBOSITY: info('Retrieving %d sequence(s)' % len(fetchme))

	clean_fetch(fetchme, outdir + '/sequences', force=force, email=email)

	if VERBOSITY: info('Done retrieving %d sequences' % len(fetchme))

	#prepare correspondences for identifind (marks B, C)
	allseqs = []
	bars = []
	seqs = {}
	pars = []
	if VERBOSITY: info('Aligning subsequences to sequences (x%d)' % len(fulltrans))
	for i, pair in enumerate(fulltrans):
		[allseqs.append(x) for x in pair]

		#bar A
		#bars.append(pairstats[pair[1]][pair[0]][3])
		#pars.append(pairstats[pair[1]][pair[0]][2])
		bars.append(pairstats[pair[1]][pair[0]][2])
		pars.append(pairstats[pair[1]][pair[0]][3])
		#bar B, C

		try: seqb = seqs[pair[1]]
		except KeyError:
			with open('%s/sequences/%s.fa' % (outdir, pair[1])) as f: seqb = seqs[pair[1]] = f.read()

		try: seqc = seqs[pair[2]]
		except KeyError:
			with open('%s/sequences/%s.fa' % (outdir, pair[2])) as f: seqc = seqs[pair[2]] = f.read()

		
		if DEBUG: info('Performing 2 subsequence-sequence alignments')
		bars.append(identifind(alnregs[pair[1]][pair[2]][0], seqb)[2:4])
		bars.append(identifind(alnregs[pair[1]][pair[2]][1], seqc)[2:4])

		#bar D
		#bars.append(pairstats[pair[2]][pair[3]][3])
		#pars.append(pairstats[pair[2]][pair[3]][2])
		bars.append(pairstats[pair[2]][pair[3]][2])
		pars.append(pairstats[pair[2]][pair[3]][3])

		try: subseqs = alnregs[pair[1]][pair[2]]
		except KeyError: subseqs = alnregs[pair[2]][pair[1]]

	#make graphs for all individual full-lengthers
	if VERBOSITY: info('Generating QUOD plots')

	for x in allseqs:
		try: seqs[x]
		except KeyError: 
			with open('%s/sequences/%s.fa' % (outdir, x)) as f: seqs[x] = f.read()

	
	for i in range(0, len(allseqs), 4):
		quod_set(tuple(allseqs[i:i+4]), seqs, outdir + '/sequences', outdir + '/graphs/', dpi=dpi, force=force, bars=bars[i:i+4], silent=not i, pars=pars[i//2:i//2+2])

	#make graphs for all pairs of sequences
	for s1 in alnregs: 
		for s2 in alnregs[s1]: 
			quod.what(alnregs[s1][s2], labels=[s1,s2], title='%s (red) vs %s (blue)' % (s1,s2), imgfmt='png', directory=outdir+'/graphs', filename='%s_vs_%s.png' % (s1,s2), dpi=dpi, hide=1)

	if VERBOSITY: info('Generating TCBLAST plots')
	blasts = {}
	tmcount = {}
	seqbank = {}
	for pair in fulltrans:
		#blasts[tuple(pair)] = [blastem(pair[1], indir=outdir, outdir=outdir, dpi=dpi), blastem(pair[2], indir=outdir, outdir=outdir, dpi=dpi, force=force, seqbank=seqbank, tmcount=tmcount, maxhits=maxhits)]
		blasts[tuple(pair)] = [blastem(pair[i+1], indir=outdir, outdir=outdir, dpi=dpi, maxhits=maxhits) for i in range(2)]

	if fulltrans:
		if VERBOSITY: info('Generating %d HTML reports' % len(fulltrans))
		for i, pair in enumerate(fulltrans):
			pairseqs = []
			for seq in pair: 
				try: pairseqs.append(seqs[seq])
				except KeyError: 
					with open('%s/sequences/%s.fa' % (outdir, seq)) as f: pairseqs.append(f.read())

			if i > 0: lastpair = fulltrans[i-1]
			else: lastpair = None
			if i < (len(fulltrans)-1): nextpair = fulltrans[i+1]
			else: nextpair = None
			build_html(pair + tuple(pairseqs) + tuple(stats[pair[1]][pair[2]]), indir=outdir, blasts=blasts[tuple(pair)], outdir=(outdir + '/html'), filename='%s_vs_%s.html' % tuple(pair[1:3]), lastpair=lastpair, nextpair=nextpair)
	else:

		if minz is None: zmin = '-inf' 
		else: zmin = '%0.1f' % minz
		if maxz is None: zmax = '+inf'
		else: zmax = '%0.1f' % maxz

		info('Generated 0 HTML reports: No significant Protocol2 hits found with Z-scores between %s and %s' % (zmin, zmax))

if __name__ == '__main__':

	import argparse

	parser = argparse.ArgumentParser(description='HTML Visualization of Reasonable, Decent Alignment Networks')

	parser.add_argument('--p1d', metavar='PATH', default=['.'], nargs='+', help='famXpander directories or table(s) (generally psiblast.tbl). Note: Running "cut -f1-12" on psiblast.tbl will greatly improve performance, but compatibility with famXpander/9.X.99/psiblast.tbl directory structures is implemented. Directory traversal is not implemented yet.')
	parser.add_argument('--p2d', metavar='PATH', default='.', help='Protocol2 directory or results table (generally results.tbl). If using on root Protocol2 directories, -f is required.')

	parser.add_argument('-o', '--outdir', metavar='DIR', default='hvordan_out', help='output directory {default:hvordan_out}')

	parser.add_argument('-f', '--fams', metavar='FAMILY', default=None, nargs=2, help='families to inspect. Required if using --p2d on root Protocol2 directories')

	parser.add_argument('-z', '--z-min', default=15, type=int, help='minimum Z score {default:15}')
	parser.add_argument('-Z', '--z-max', default=None, type=int, help='maximum Z score {default:none}')

	parser.add_argument('-c', '--clobber', action='store_true', help='force redownloads/regenerates where applicable')
	parser.add_argument('-r', '--dpi', type=int, default=100, help='resolution of graphs {default:100}')
	parser.add_argument('-m', '--max-hits', type=int, default=10, help='how many TCBLAST hits to BLAST for. Contributes significantly to execution time for small famXpander results. {default:10}')

	if 'ENTREZ_EMAIL' in os.environ:
		parser.add_argument('-e', '--email', default=None, help='Working email in case too many requests get sent and the NCBI needs to initiate contact. Defaults to checking $ENTREZ_EMAIL if set. {current value: %s}' % os.environ['ENTREZ_EMAIL'])
	else: parser.add_argument('-e', '--email', default=None, help='Working email in case too many requests get sent and the NCBI needs to initiate contact. Defaults to checking $ENTREZ_EMAIL if set. {unset}')

	parser.add_argument('-i', metavar='ACC', nargs='+', help='Operate only on pairs containing these accessions')
	parser.add_argument('-p', metavar='ACC', nargs='+', help='Operate only on these specific pairs.')

	args = parser.parse_args()

	if args.p1d == '.' and args.p2d == '.': 
		parser.print_help()
		exit()

	summarize(args.p1d, args.p2d, args.outdir, minz=args.z_min, maxz=args.z_max, dpi=args.dpi, force=args.clobber, email=args.email, musthave=args.i, thispair=args.p, fams=args.fams, maxhits=args.max_hits)
