import numpy as np
import pandas as pd
import os
import regex
import Bio
import Bio.SeqIO
from bisect import bisect_left
import pickle
from CAI import CAI
from Bio.Blast import NCBIWWW, NCBIXML


class Match:
	def __init__(self, string, start, stop):
		self.string = string
		self.start_pos = start
		self.stop_pos = stop
		self._group0 = string[start:stop]
	def group(self):
		return self._group0
	def span(self):
		return (self.start_pos, self.stop_pos)
	def __repr__(self):
		return f"<Match object; span={self.span()}, match='{self.group()}'>"

def orf_iter(sequence):
	# Cerco tutti i codoni di stop e salvo il punto di partenza
	stops = []
	for stop in regex.finditer(r"TAA|TAG|TGA", sequence):
		stops.append(stop.span()[0])
	# Itero sugli start
	for start in regex.finditer(r"ATG", sequence):
		# cerco il primo stop con valore di distanza multiplo di 3
		pos_start = start.span()[0]
		pos_stop = bisect_left(stops, pos_start + 3)
		while (pos_stop < len(stops)) and ((stops[pos_stop] - pos_start) % 3 != 0):
			pos_stop = pos_stop + 1
		if pos_stop < len(stops):
			yield Match(sequence, pos_start, stops[pos_stop])
			
################################################################################
#                         ALGORITMO DI GENE PREDICTION                         #
# Scopo del progetto è, data la sequenza scaricata, salvare in un dataframe    #
# la posizione delle ORF presenti nella sequenza e l'indicazione se si tratta  #
# di regione codificante o non codificante e se è stata verificata su BLAST    #
################################################################################

#caricamento sequenza
seq = next(Bio.SeqIO.parse(f"group_31_seq.fasta", "fasta")).seq
#setup df
results = pd.DataFrame(columns=["ORFpos","Coding","BLAST"])

# Calcolo frequenze CG in tutta la sequenza
exp_cpg = seq.count("C")/len(seq) * seq.count("G")/len(seq)

#funzioni per sensori di segnale
def tatabox(seq):
    if regex.search(r"TATA[AT]A[AT]", str(seq)):
        return True
    return False
    
def inr(seq):
    if regex.search(r"[CT][CT]A[AGCT][AT][CT][CT]", str(seq)):
        return True
    return False
      
def cpg(seq):
    seq_monte = str(seq[orf.span()[0]-3000:orf.span()[0]])
    #contenuto G e G
    GC_monte = len(regex.findall(r"[GC]", seq_monte))/3000
    #per prima cosa la frequenza deve essere superiore a 0.5
    if GC_monte > 0.5:
        #rapporto osservato/atteso
            obs_CG = seq_monte.count("CG")/3000
            #inoltre il rapporto osservato/atteso deve essere maggiore di 0.6
            if obs_CG/exp_cpg > 0.6:
                #se è così aggiunto 1 punto
                return True
    return False

def kozak(seq):
    if regex.search(r"[AG]CCATGG", str(seq)):
        return True
    return False

#funzioni per sensori di contenuto
def gc_content(seq):
    GC_content = len(regex.findall(r"[GC]", str(seq)))/3000
    if GC_content > 0.45 and GC_content < 0.65:
        return True
    return False

def codon_freq(seq):
    cai_dict = {"GCT": 1.000, "TTG": 0.600, "GCC": 0.900, "CTA": 0.400, "GCA": 0.700, "TTA": 0.300, "GCG": 0.600, "AAG": 1.000, "CGT": 1.000, "AAA": 0.800, "CGC": 0.900, "ATG": 1.000, "AGA": 0.700, "TTC": 1.000, "CGG": 0.600, "TTT": 0.800, "AGG": 0.500, "CCG": 1.000, "CGA": 0.400, "CCC": 0.900, "AAC": 1.000, "CCT": 0.800, "AAT": 0.800, "CCA": 0.700, "GAC": 1.000, "AGC": 1.000, "GAT": 0.800, "TCT": 0.900, "TGC": 1.000, "TCC": 0.800, "TGT": 0.800, "AGT": 0.700, "CAG": 1.000, "TCA": 0.600, "CAA": 0.700, "TCG": 0.500, "GAG": 1.000, "ACC": 1.000, "GAA": 0.800, "ACT": 0.900, "GGC": 1.000, "ACA": 0.700, "GGT": 0.900, "ACG": 0.600, "GGA": 0.700, "TGG": 1.000, "GGG": 0.600, "TAC": 1.000, "CAC": 1.000, "TAT": 0.800, "CAT": 0.800, "GTG": 1.000, "ATC": 1.000, "GTC": 0.900, "ATT": 0.900, "GTT": 0.800, "ATA": 0.500, "GTA": 0.600, "CTG": 1.000, "TAA": 1.000, "CTC": 0.800, "TGA": 0.900, "CTT": 0.700, "TAG": 0.800}
    cai_score=CAI(seq, weights=cai_dict)
    if cai_score > 0.8:
        return True
    return False

def verifica_sequenza_blast(seq):
    try:
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="refseq_rna",
            sequence=seq
        )

        blast_records = NCBIXML.read(result_handle)
        
        # Controlla se ci sono match con E-value < 1e-6
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 1e-6:
                    return True
        return False
    except Exception as e:
        print(f"Errore durante l'esecuzione di BLAST: {e}")
        return False

#VARIABILI SCORE
score_tatabox = 3
score_inr = 2
score_cpg = 1
score_kozak = 2
score_gc = 1
score_codon = 1

#VARIABILI POSIZIONALI
inizio_tatabox = -40
fine_tatabox = -20
inizio_inr = -10
fine_inr = 10
inizio_cpg = -3000
inizio_kozak = -10
fine_kozak = 10


#inizio analisi
for orf in orf_iter(str(seq)):
  # FILTRAGGIO PRELIMINARE
    if orf.span()[1] - orf.span()[0] < 150:
        continue
    results.loc[len(results)] = [orf.span()[0], False, False]

    #PARTE CODING
    #inizializzo score
    score = 0
    #il tatabox si trova tra -35 e -25 prima dell'ATG iniziale
    subseq = str(seq[orf.span()[0]+inizio_tatabox:orf.span()[0]+fine_tatabox])
    if tatabox(subseq):
        score = score + score_tatabox
    subseq = str(seq[orf.span()[0]+inizio_inr:orf.span()[0]+fine_inr])
    #l'inr si trova tra -10 e +10 attorno all'ATG iniziale
    if inr(subseq):
        score = score + score_inr
    #solo se la sequenza inizia dopo 3000 basi
    if orf.span()[0] > 3000:
        subseq = str(seq[orf.span()[0]+inizio_cpg:orf.span()[0]])
        #controllo se nelle 3000 basi precedenti c'è un'alta frequenza di CpG
        if cpg(subseq):
            score = score + score_cpg
    #il kozak si trova tra -10 e +10 attorno all'ATG iniziale
    subseq = str(seq[orf.span()[0]+inizio_kozak:orf.span()[0]+fine_kozak])
    if kozak(subseq):
        score = score + score_kozak
    #i sensori di contenuto si applicano alla sequenza da ATG a stop
    subseq = str(seq[orf.span()[0]:orf.span()[1]])
    if gc_content(subseq):
        score = score + score_gc
    subseq = str(seq[orf.span()[0]:orf.span()[1]])
    if codon_freq(subseq):
        score = score + score_codon
    if score >= 5:
        results.loc[len(results)-1, 'Coding'] = True
        #PARTE VERIFICA BLAST
        print(f'Controllo se la sequenza in posizione {orf.span()} è verificata su BLAST')
        results.loc[len(results)-1, 'BLAST'] = verifica_sequenza_blast(subseq)

#filtra solo le ORF codificanti
cod = results[results['Coding'] == True]
print(cod)

#filtra solo le ORF codificanti verificate su BLAST
ver = results[results['BLAST'] == True]
print(ver)
