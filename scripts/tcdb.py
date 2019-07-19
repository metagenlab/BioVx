#!/usr/bin/env python

# A TCDB Interface
import urllib
import tempfile
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
#### the program does not seeem to use the mysql interface at all,
#### so maybe should be out for now?
#### I agree, all database functions are removed -V
from Bio import SeqIO
from datetime import datetime
import pickle

def acc2fasta(accs):
    url='http://www.tcdb.org/projectv/acc2fasta.php?accs=%s'%(','.join(accs))
    result = urllib.urlopen(url)
    text = result.read()
    fastafile=tempfile.TemporaryFile()
    fastafile.write(text)
    fastafile.seek(0)
    myfastas = tempfile.NamedTemporaryFile()
    fastas = SeqIO.parse(fastafile,'fasta')
    for fasta in fastas:
        head = fasta.description.split("@@@")
        acc = head[1]
        fasta.id=acc
        SeqIO.write(fasta,myfastas,'fasta')
    myfastas.seek(0)
    return myfastas

def define_family(family,fasta=False):
    url = "http://www.tcdb.org/projectv/sft.php?fam=%s" %(family)
    response = urllib.urlopen(url)
    acc = response.read().strip()
    accs = [i.strip() for i in acc.split(' ')]
    if accs[0] == '':
        return False
    if fasta is False:
        return accs
    # Get Seq Object instead
    accs = ",".join(accs)
    url = 'http://www.tcdb.org/cgi-bin/projectv/accs2fasta.py?accs=%s'%accs
    response = urllib.urlopen(url)
    fastas = SeqIO.parse(response,'fasta')
    return fastas

def download_fasta(out):
    urllib.urlretrieve ('http://www.tcdb.org/api.php?tcid=all', out)
    return

def prep_db(path):
    locald = '/'.join(path.split('/')[0:-1])
    if not os.path.exists(locald):
        os.makedirs(locald)
    urllib.urlretrieve("http://tcdb.org/public/tcdb", path)
    os.system('makeblastdb -dbtype prot -in '+path)

def prep_bb(path):
    urllib.urlretrieve("http://tcdb.org/public/betabarrel", path)

def use_local(path, update=True):
    if os.path.exists(path) is False:
        prep_db(path)
        return
    if update:
        stat = os.stat(path)
        fileage = datetime.fromtimestamp(stat.st_mtime)
        now = datetime.now()
        delta = now-fileage
        if delta.days >= 5:
            os.remove(path)
            prep_db(path)
            return
    return

def use_local_betabarrel(path, update=True):
    if os.path.exists(path) is False:
        prep_bb(path)
        return
    if update:
        stat = os.stat(path)
        fileage = datetime.fromtimestamp(stat.st_mtime)
        now = datetime.now()
        delta = now-fileage
        if delta.days >= 5:
            os.remove(path)
            prep_db(path)
            return
    return

class Names:
    # This wrapper is for handling all nomenclature.
    def __init__(self):
        # Load Abbreviations for families
        self.familyabr = {}
        abr = urllib.urlopen('http://tcdb.org/cgi-bin/projectv/family_abbreviations.py').read().split("\n")
        for a in abr:
            try:
                (family,symbol) = a.split('\t')
                symbol = symbol if bool(len(symbol)) else family
                self.familyabr.setdefault(family,symbol)
            except:
                pass

    def get_family_abr(self,family):
        if family in self.familyabr.keys():
            return self.familyabr[family]
        return family

class Substrates:
    # This wrapper allows access to TCID->Substrate relations
    def __init__(self):
        #data=urllib.urlopen('http://tcdb.org/cgi-bin/projectv/substrates.py').read().split('\n')
        self.mysubstrates = {}
        #for line in data[:-1]:
        #(tcid,substrate) = line.split('\t')
        #self.mysubstrates.setdefault(tcid,[]).append(substrate)

    def get_tcid_substrates(self,tcid):
        
        substrates = pickle.load(open(os.environ['PYTHONPATH']+'/substrate_dict.py','rb'))

        if tcid in substrates:
		
            return substrates[tcid][1:6]

        else:

            return None





