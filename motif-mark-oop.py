#!/usr/bin/env python

# allow me to use Exon in a type annotation before define class Exon 
from __future__ import annotations

import argparse
import re
import cairo

def get_args():
    parser = argparse.ArgumentParser(description="This code will parse through a fasta file and a motif file and it will return an PNG image of the intron, exon and motif locations.")
    parser.add_argument("-f",  "--fasta", help="Absolute file path to FASTA file.", required = True)
    parser.add_argument("-m", "--motifs", help = "Absolute file path to the motifs file (max 5 motifs).", required = True)
    return parser.parse_args()

args = get_args()
fasta_file = args.fasta
motif_file = args.motifs

def oneline_fasta(file: str, new_file: str):
    '''This function takes a fasta file with multiple sequence lines and combine them together.'''
    seq=""
    header = ""
    new_fh = open(new_file, "w")
    with open(file, "r") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if seq != "":
                    print(f"{header}\n{seq}", file=new_fh)
                header = line
                seq = "" ##restart sequence at each new record
            else:
                seq+=line
    print(f"{header}\n{seq}", file=new_fh)
    new_fh.close()
    return None

###########
# Classes #
###########

class Gene:
    '''This is how you create a gene object that draw itself.'''
    def __init__(self, seq, y, x_shift = 30): ## setting default y to 30 so won't start drawing at origin
        ## Data ##
        self.x_shift = x_shift 
        self.seq = seq ## saving sequence
        self.x_start = 0## starting at origin
        self.x_end = len(seq) ## length of the sequence is the x end
        self.y = y # y1 and y2 are the same because we are drawing a line
        self.exon = self.create_exon() ## calling the function to create an exon
        #self.text = text ## text to write

    def __repr__(self):
        '''Changes the representation of the object.''' 
        return f"Gene('{self.seq}', {self.y})"

    ## Methods (functions) ## 

    def create_exon(self) -> Exon:
        '''This function parse through the sequence line in a fasta file and 
        finds the exon and its coordinates.'''
        match = re.search("[ATCGU]+", self.seq) ## it is a Match (class)
        if match is not None:
            x_start = match.start() ## start of the match
            x_end = match.end() - 1 ## end of the match (removing 1)
        return Exon(x_start, x_end, self.y)


    def draw(self, ctx):
        ''' This function draw a line with the size of the gene.'''
        ctx.set_line_width(2) ## intron line
        ctx.set_source_rgb(0, 0, 0) ## setting color to black
        ctx.move_to(self.x_start + self.x_shift, self.y)
        ctx.line_to(self.x_end + self.x_shift, self.y) 
        ctx.stroke()
        print("Line draw!")
        #print(type(ctx))
    
    def write_header(self, ctx, text):
        ''' This function include the header for each gene. '''
        match = re.search(">(.+) chr(.+):(.+)-(\d+)", text)
        if match is not None:
            gene_name = match.group(1) ## first match is the gene
            chr = match.group(2) ## second match is chromossome
            loc_start = match.group(3) ## third match is start location
            loc_end = match.group(4) ## fourth match is end location
            ctx.set_source_rgb(0, 0, 0) 
            ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            ctx.set_font_size(15)
            ctx.move_to(self.x_start + 10, self.y - 35)
            ctx.show_text(f"Gene: {gene_name}    Chromossome: {chr}    Location: {loc_start}-{loc_end}")
            ctx.stroke()

class Exon:
    ''' This is how you create a exon object. ''' 
    def __init__(self, x_start:int, x_end:int, y:int, x_shift = 30):
        self.x_start = x_start
        self.x_end = x_end
        self.y = y
        self.x_shift = x_shift 

    def __repr__(self):
        '''Changes the representation of the object.''' 
        return f"Exon({self.x_start}, {self.x_end}, {self.y})"
    
    ## Methods (functions) ## 

    def draw(self, ctx):
        ''' This function draw an exon using their start and end position. '''
        ctx.set_line_width(40)
        ctx.set_source_rgb(0, 0, 0) ## setting color to black
        ctx.move_to(self.x_start+self.x_shift, self.y) ## adding the start position of the line
        ctx.line_to(self.x_end+self.x_shift, self.y) ## adding the start position of the line
        ctx.stroke()
        print("Exon draw!")


class Motif:
    ''' This is how you create a motif object.'''
    def __init__(self, seq, motif_string, y, x_shift = 30): ## setting default y to 30 so won't start drawing at origin
        ## Data ##
        self.seq = seq
        self.motif_string = motif_string ## motifs before convertion (e.g. YYYY)
        self.amb_regex = self.get_ambiguous_regex() ## calling the function to convert nucleotides
        self.y = y
        self.x_shift=x_shift
        #self.color = color
        #self.mot_find = self.find_motifs(seq) 

    ## Methods ##
    def get_ambiguous_regex(self):
        '''This function parse through the sequence line in a fasta file and 
        convert ambiguous nucleotides into their IUPAC equivalents.'''

        ## creating a dictionary to store degenerate nucleotides possibilities (key: nucletide, value: list with bases)
        iupac_nucleotides_dict = {"A":"[Aa]", 
        "C": "[Cc]", 
        "G": "[Gg]", 
        "T": "[Tt]", 
        "U": "[Uu]", 
        "W": "[AaTt]", ##weak
        "S": "[CcGg]", ##strong
        "M": "[AaCc]", ##amino
        "K": "[GgTt]", ##ketone
        "R": "[AaGg]", ##purine
        "Y": "[CcTt]", ##pyrimidine
        "B": "[CcGgTt]", ##not A
        "D": "[AaGgTt]", ##not C
        "H": "[AaCcTt]", ##not G
        "V": "[AaCcGg]", ##not T
        "N": "[AaCcGgTt]", ##any one base
        }

        regex = "" ## creating empty string to hold the regex sequence
        motif = self.motif_string.upper() ## convert in upper case

        for letter in motif:
            if letter in iupac_nucleotides_dict:
                regex += iupac_nucleotides_dict[letter] ## letter is now the regex
            else:
                regex += letter ## without conversion
        return regex

    def find_motifs(self, seq:str):
        '''This function parse through a fasta sequence and finds all motifs that match the motif list after 
        ambiguous transformation.'''
        regex = self.amb_regex
        motif_match = {}
        #sten_list = [] ## list to hold start and end positions
        for m in re.finditer(rf"{regex}", seq.upper()): ## r for regex and f for f-string
            if m is not None:
                if m.group() not in motif_match:
                    #print(m.group(), m.start())
                    motif_match[m.group()] = [(m.start(), m.end())]
                else:
                    motif_match[m.group()].append((m.start(), m.end())) ## appending list with more than one occurrence of the motif
                #motif_match[m.group()]=(m.start(), m.end()) ## saving match motif in a dictionary with the start and end positions
        return motif_match
    
    def draw(self, ctx, rgb): ## color as a positional argument
        ''' This function draws a motif. Each motif has a different color.'''
        motifs = self.find_motifs(self.seq) ## call find_motifs function inside the draw to have the start and end positions
        print(motifs)
        for motif_name, list_of_tuples in motifs.items():
            #print(motif_name, list_of_tuples)
            for tup in list_of_tuples:
                x_start = tup[0]
                x_end = tup[1]
                ctx.set_line_width(40) ## the same as the exon
                ctx.set_source_rgb(rgb[0], rgb[1], rgb[2]) 
                ctx.move_to(x_start + self.x_shift, self.y) ## adding the start position of the line
                ctx.line_to(x_end + self.x_shift, self.y) ## adding the start position of the line
                ctx.stroke()
                ctx.set_source_rgb(0,0,0)
            print("Motifs draw!")


########
# Main #
########

oneline_fasta("Figure_1.fasta", "Figure1_one_line.fasta") ## output file with one line sequence

## creating motif list
motif_list = []
with open(motif_file, "r") as mf:
    for line in mf:
        line = line.strip("\n").upper()
        if "U" in line:
            line = line.replace("U", "T") 
            motif_list.append(line)
        else:
            motif_list.append(line)
    print(motif_list)

## loop to find the height of the canvas
lengths = []
g = 0 ## counts for the number of genes
with open(fasta_file, "r") as ff:
    for line in ff:
        if not line.startswith(">"):
            lengths.append(len(line))
        else:
            g +=1

## setting canvas 
height = g * 170  ## depending on the number of genes in the file
width = max(lengths) + 200 ## depending on the size of the sequences in the fasta
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height) ## setting surface to png
ctx = cairo.Context(surface) ## setting context
ctx.set_source_rgb(1, 1, 1) ## white background
ctx.paint()

with open(fasta_file, "r") as fq:  ## adjust this to move to the other line and draw there
    i = 0
    color_list = [(0.39, 0.58, 0.92), (0.85, 0.43, 0.83), (0.12, 0.69, 0.66), (0.57,0.43,0.85), (1, 0.84, 0)] ## blue, pink, green, purple, yellow
    while True:
        #y = 10 ## increase y for each gene !! 
        header = fq.readline()
        seq = fq.readline()
        i +=1 ## counter for the lines
        if header == "": 
            break   ## EOF (end of file)
        gene = Gene(seq, i*100)
        gene.write_header(ctx, header)
        gene.draw(ctx)
        #print(gene)
        exon = gene.create_exon()
        exon.draw(ctx)
        for j, mname in enumerate(motif_list): ## motif index (to match the color) and motif name
            motif = Motif(seq, mname, i*100) ## changing y with every line
            motif.draw(ctx, color_list[j]) 
l_y_start = i*100
l_x_start = 30

## including legend
ctx.set_source_rgb(0,0,0)
ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
ctx.move_to(l_x_start, l_y_start - 20)
ctx.set_font_size(20)
ctx.show_text("Motifs")
ctx.set_font_size(12)
for j, mname in enumerate(motif_list): ## draw legend for each motif
    ctx.set_source_rgb(color_list[j][0], color_list[j][1], color_list[j][2]) ## setting color 
    ctx.rectangle(l_x_start, l_y_start, 20, 20) ## (x0, y0, w, h)
    ctx.fill()
    ctx.set_source_rgb(0,0,0)
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.move_to(l_x_start + 25, l_y_start + 15)
    ctx.show_text(mname)  ## motif name
    ctx.stroke()
    l_y_start += 30
surface.write_to_png(str(fasta_file) + ".png") ## saving file




